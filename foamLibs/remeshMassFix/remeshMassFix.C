/*
Copyright (c) 2026 Zhi Sheng Chen

Licensed under the PolyForm Noncommercial License 1.0.0 (the "License").
You may not use this file except in compliance with the License. Any
commercial use requires a separate commercial license from the copyright
holder.

You may obtain a copy of the License at:
    https://polyformproject.org/licenses/noncommercial/1.0.0

Or see the LICENSE file in the root of this repository.

This software is provided "as is", without warranty of any kind, express or
implied. See the License for the specific language governing permissions
and limitations.
*/

// remeshing helper: map fix and mass fix after mapping onto a new mesh
//   remeshMassFix -time <t>                prints "gasMass <m>" (when rho exists)
//   remeshMassFix -time <t> -writeCover    writes remeshCover = 1 on the old mesh before mapping, and fills the
//                                          matrix temperature Ts_<zone> of each porous zone outside the zone with
//                                          the zone values (a few cell layers): outside the zone Ts stays at its
//                                          start value, the map would mix it into the zone's edge cells
//   remeshMassFix -time <t> -combine       after the volume-weighted map (mapFieldsPar, cellVolumeWeight, weights
//                                          normalised over the overlap): new cells with no overlap at all are left
//                                          at zero (snapped walls never match exactly, T = 0 crashes the thermo).
//                                          They are found from remeshCover (< -minCover, default 0.5) and filled
//                                          from their covered neighbours, layer by layer (across processor
//                                          patches); processor patch values are then evaluated from the cells.
//                                          The fill replaces a second, nearest-cell map: v2412 mapNearest crashes
//                                          in parallel (cell walk over a disconnected piece of the old mesh)
//   remeshMassFix -time <t> -mass <m>      scales p and rho by m/gasMass (T kept, rho ~ p for a perfect gas)
// setExprFields would do the scaling, but v2412 aborts on any expression that reads a field

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "IOobjectList.H"
#include "calculatedFvPatchFields.H"
#include "OSspecific.H"
#include "IOmanip.H"
#include "syncTools.H"

using namespace Foam;

static const word coverName("remeshCover");


// cells not set yet = average of their set neighbours, one layer per pass (across processor patches too), until
// every cell is set or maxLayers passes are done. returns the number of cells filled
template<class Type>
label fillFromNeighbours
(
    const fvMesh& mesh,
    Field<Type>& values,
    scalarField& isSet,
    const label maxLayers
)
{
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const label nInternal = mesh.nInternalFaces();
    label nFilled = 0;

    for (label layer = 0; layer < maxLayers; ++layer)
    {
        Field<Type> sum(mesh.nCells(), Zero);
        scalarField count(mesh.nCells(), 0);
        for (label facei = 0; facei < nInternal; ++facei)
        {
            const label o = own[facei];
            const label n = nei[facei];
            if (isSet[n] > 0.5 && isSet[o] < 0.5) { sum[o] += values[n]; count[o] += 1; }
            if (isSet[o] > 0.5 && isSet[n] < 0.5) { sum[n] += values[o]; count[n] += 1; }
        }
        // neighbours on the other side of processor patches
        scalarField nbrSet(mesh.nBoundaryFaces());
        Field<Type> nbrValue(mesh.nBoundaryFaces());
        syncTools::swapBoundaryCellList(mesh, isSet, nbrSet);
        syncTools::swapBoundaryCellList(mesh, values, nbrValue);
        for (const polyPatch& pp : mesh.boundaryMesh())
        {
            if (!pp.coupled())
            {
                continue;
            }
            forAll(pp, i)
            {
                const label bFacei = pp.start() + i - nInternal;
                const label o = own[pp.start() + i];
                if (nbrSet[bFacei] > 0.5 && isSet[o] < 0.5) { sum[o] += nbrValue[bFacei]; count[o] += 1; }
            }
        }
        label nLayer = 0;
        forAll(values, celli)
        {
            if (count[celli] > 0)
            {
                values[celli] = sum[celli]/count[celli];
                isSet[celli] = 1;
                ++nLayer;
            }
        }
        nFilled += nLayer;
        if (returnReduce(nLayer, sumOp<label>()) == 0)
        {
            break;
        }
    }
    return nFilled;
}


template<class Type>
label combineFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    const volScalarField& cover,
    const scalar minCover,
    label& nUncovered,
    label& nUnfilled
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    label nFields = 0;

    for (const word& name : objects.sortedNames(fieldType::typeName))
    {
        if (name == coverName)
        {
            continue;
        }
        fieldType fld
        (
            IOobject(name, mesh.time().timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
            mesh
        );

        // covered cells: volume-weighted average over the overlap, the others from their neighbours
        Field<Type>& values = fld.primitiveFieldRef();
        scalarField isSet(mesh.nCells(), 0);
        label nLocal = 0;
        forAll(values, celli)
        {
            if (cover[celli] >= minCover)
            {
                values[celli] /= cover[celli];
                isSet[celli] = 1;
            }
            else
            {
                ++nLocal;
            }
        }
        fillFromNeighbours(mesh, values, isSet, 20);
        label nLeft = 0;
        forAll(isSet, celli)
        {
            if (isSet[celli] < 0.5)
            {
                ++nLeft;
            }
        }
        nUncovered = returnReduce(nLocal, sumOp<label>());
        nUnfilled = returnReduce(nLeft, sumOp<label>());

        // processor (coupled) patches from the filled cells: the mapped values next to uncovered cells are zero,
        // and the thermo is built from them before the solver updates anything (T = 0 -> FPE). Wall values come
        // from the patch map (nearest face)
        fld.boundaryFieldRef().evaluate_if
        (
            [](const fvPatchField<Type>& pfld) { return pfld.coupled(); }
        );

        fld.write();
        ++nFields;
    }
    return nFields;
}


// porous zones with a matrix temperature: Ts_<zone> fields whose cellZone exists
wordList matrixFields(const fvMesh& mesh, const IOobjectList& objects)
{
    wordList names;
    for (const word& name : objects.sortedNames(volScalarField::typeName))
    {
        if (name.size() > 3 && name.substr(0, 3) == "Ts_" && mesh.cellZones().findZoneID(name.substr(3)) != -1)
        {
            names.append(name);
        }
    }
    return names;
}


// zone volume and volume-averaged matrix temperature
void matrixReport(const fvMesh& mesh, const volScalarField& Ts, const word& when)
{
    const labelList& cells = mesh.cellZones()[mesh.cellZones().findZoneID(Ts.name().substr(3))];
    scalar vol = 0;
    scalar sumTV = 0;
    for (const label celli : cells)
    {
        vol += mesh.V()[celli];
        sumTV += Ts[celli]*mesh.V()[celli];
    }
    reduce(vol, sumOp<scalar>());
    reduce(sumTV, sumOp<scalar>());
    Info<< "matrix " << Ts.name() << " " << when << ": zone volume " << setprecision(10) << vol
        << " m3, mean Ts " << sumTV/max(vol, VSMALL) << " K" << endl;
}


// Ts outside its zone from the zone values, nLayers cell layers deep
void extendMatrixField(const fvMesh& mesh, volScalarField& Ts, const label nLayers)
{
    scalarField isSet(mesh.nCells(), 0);
    for (const label celli : mesh.cellZones()[mesh.cellZones().findZoneID(Ts.name().substr(3))])
    {
        isSet[celli] = 1;
    }
    const label nFilled = returnReduce(fillFromNeighbours(mesh, Ts.primitiveFieldRef(), isSet, nLayers), sumOp<label>());
    Info<< "matrix " << Ts.name() << ": " << nFilled << " cells around the zone set from the zone (" << nLayers
        << " layers)" << endl;
}


int main(int argc, char *argv[])
{
    argList::addNote("remeshing: gas mass, coverage field, combined map, mass fix");
    argList::addOption("mass", "kg", "gas mass to restore: p and rho are scaled by mass/gasMass");
    argList::addBoolOption("writeCover", "write remeshCover = 1 (old mesh, before mapping)");
    argList::addBoolOption("combine", "fill the cells the volume-weighted map left without overlap");
    argList::addOption("minCover", "frac", "below this covered fraction a cell is filled from its neighbours (0.5)");
    timeSelector::addOptions(false, false);

    #include "setRootCase.H"
    #include "createTime.H"

    const instantList times = timeSelector::select0(runTime, args);
    if (times.size() != 1)
    {
        FatalErrorInFunction
            << "select one time (-time <t> or -latestTime)" << exit(FatalError);
    }
    runTime.setTime(times[0], 0);

    // mesh of that time (a remeshed time directory holds its own polyMesh)
    fvMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    if (args.found("writeCover"))
    {
        volScalarField cover
        (
            IOobject(coverName, runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh,
            dimensionedScalar(dimless, 1),
            calculatedFvPatchScalarField::typeName
        );
        cover.write();
        Info<< "wrote " << coverName << endl;

        const IOobjectList objects(mesh, runTime.timeName());
        for (const word& name : matrixFields(mesh, objects))
        {
            volScalarField Ts
            (
                IOobject(name, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
                mesh
            );
            matrixReport(mesh, Ts, "before the map");
            extendMatrixField(mesh, Ts, 4);
            Ts.write();
        }
    }

    if (args.found("combine"))
    {
        const scalar minCover = args.getOrDefault<scalar>("minCover", 0.5);
        const IOobjectList objects(mesh, runTime.timeName());
        label nFields = 0;
        label nUncovered = 0;
        label nUnfilled = 0;
        {
            const volScalarField cover
            (
                IOobject(coverName, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
                mesh
            );
            nFields += combineFields<scalar>(mesh, objects, cover, minCover, nUncovered, nUnfilled);
            nFields += combineFields<vector>(mesh, objects, cover, minCover, nUncovered, nUnfilled);
            Info<< "coverage of the new cells by the old mesh: min " << gMin(cover.primitiveField())
                << ", max " << gMax(cover.primitiveField()) << endl;
            rm(cover.objectPath());
        }
        for (const word& name : matrixFields(mesh, IOobjectList(mesh, runTime.timeName())))
        {
            const volScalarField Ts
            (
                IOobject(name, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
                mesh
            );
            matrixReport(mesh, Ts, "after the map");
        }
        Info<< "combined " << nFields << " fields, " << nUncovered << " cells without overlap filled from their "
            << "neighbours (covered below " << minCover << ")" << endl;
        if (nUnfilled)
        {
            FatalErrorInFunction
                << nUnfilled << " cells without overlap have no filled neighbour within 20 layers" << exit(FatalError);
        }
    }

    // gas mass (compressible only)
    IOobject rhoIO("rho", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE);
    if (!rhoIO.typeHeaderOk<volScalarField>(true))
    {
        if (args.found("mass"))
        {
            FatalErrorInFunction << "no rho in " << runTime.timeName() << exit(FatalError);
        }
        Info<< "End" << endl;
        return 0;
    }
    volScalarField rho(rhoIO, mesh);
    const scalar mass0 = gSum(rho.primitiveField()*mesh.V().field());
    Info<< "gasMass " << setprecision(15) << mass0 << endl;

    if (args.found("mass"))
    {
        const scalar mass = args.get<scalar>("mass");
        const dimensionedScalar factor(dimless, mass/mass0);
        volScalarField p
        (
            IOobject("p", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
            mesh
        );
        rho *= factor;
        p *= factor;
        rho.write();
        p.write();
        Info<< "scaled p and rho by " << setprecision(15) << factor.value()
            << " (gas mass " << mass0 << " -> " << mass << ")" << endl;
    }

    Info<< "End" << endl;
    return 0;
}
