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
//   remeshMassFix -time <t> -writeCover    writes remeshCover = 1 on the old mesh before mapping
//   remeshMassFix -time <t> -combine       after the two maps (volume weighted, nearest as <field>_nearest):
//                                          field = mapped/remeshCover where the old mesh covered more than
//                                          -minCover (default 0.5) of the cell, else the nearest value; wall
//                                          values from the nearest map. The volume-weighted map (weights
//                                          normalised over the overlap) leaves new cells with no overlap at all
//                                          at zero (snapped walls never match exactly): T = 0 crashed the thermo
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

using namespace Foam;

static const word coverName("remeshCover");
static const word nearSuffix("_nearest");


bool isNearest(const word& name)
{
    return name.size() > nearSuffix.size() && name.substr(name.size() - nearSuffix.size()) == nearSuffix;
}


template<class Type>
label combineFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    const volScalarField& cover,
    const scalar minCover,
    label& nNearest
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;
    label nFields = 0;

    for (const word& name : objects.sortedNames(fieldType::typeName))
    {
        if (isNearest(name) || name == coverName || !objects.found(name + nearSuffix))
        {
            continue;
        }
        const word timeName(mesh.time().timeName());
        fieldType fld
        (
            IOobject(name, timeName, mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
            mesh
        );
        const fieldType nearFld
        (
            IOobject(name + nearSuffix, timeName, mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
            mesh
        );

        // cells: volume-weighted average over the covered part, nearest where the old mesh hardly reaches
        Field<Type>& values = fld.primitiveFieldRef();
        forAll(values, celli)
        {
            const scalar c = cover[celli];
            if (c > minCover)
            {
                values[celli] /= c;
            }
            else
            {
                values[celli] = nearFld[celli];
                ++nNearest;
            }
        }

        // walls: nearest face (the area-weighted patch map has the same coverage problem)
        forAll(fld.boundaryField(), patchi)
        {
            if (!fld.boundaryField()[patchi].coupled())
            {
                fld.boundaryFieldRef()[patchi] == nearFld.boundaryField()[patchi];
            }
        }
        // processor (coupled) patches from the combined cells: the mapped values next to uncovered cells are zero,
        // and the thermo is built from them before the solver updates anything (T = 0 -> FPE)
        fld.boundaryFieldRef().evaluate_if
        (
            [](const fvPatchField<Type>& pfld) { return pfld.coupled(); }
        );

        fld.write();
        rm(nearFld.objectPath());
        ++nFields;
    }
    return nFields;
}


int main(int argc, char *argv[])
{
    argList::addNote("remeshing: gas mass, coverage field, combined map, mass fix");
    argList::addOption("mass", "kg", "gas mass to restore: p and rho are scaled by mass/gasMass");
    argList::addBoolOption("writeCover", "write remeshCover = 1 (old mesh, before mapping)");
    argList::addBoolOption("combine", "combine the volume-weighted and the nearest map");
    argList::addOption("minCover", "frac", "below this covered fraction a cell takes the nearest value (0.5)");
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
    }

    if (args.found("combine"))
    {
        const scalar minCover = args.getOrDefault<scalar>("minCover", 0.5);
        const IOobjectList objects(mesh, runTime.timeName());
        label nNearest = 0;
        label nFields = 0;
        {
            const volScalarField cover
            (
                IOobject(coverName, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
                mesh
            );
            nFields += combineFields<scalar>(mesh, objects, cover, minCover, nNearest);
            nFields += combineFields<vector>(mesh, objects, cover, minCover, nNearest);
            Info<< "coverage of the new cells by the old mesh: min " << gMin(cover.primitiveField())
                << ", max " << gMax(cover.primitiveField()) << endl;
            rm(cover.objectPath());
        }
        // leftovers of the nearest map (the coverage field itself, fields without a volume-weighted twin)
        for (const word& name : objects.sortedNames())
        {
            if (isNearest(name))
            {
                rm(objects[name]->objectPath());
            }
        }
        reduce(nNearest, sumOp<label>());
        Info<< "combined " << nFields << " fields, " << nNearest/max(nFields, label(1))
            << " cells took the nearest value (covered below " << minCover << ")" << endl;
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
