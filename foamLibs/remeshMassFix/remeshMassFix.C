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

// gas mass of a time directory, and the mass fix after mapping onto a new mesh (remeshing):
//   remeshMassFix -time <t>                prints "gasMass <m>"
//   remeshMassFix -time <t> -mass <m>      scales p and rho by m/gasMass (T kept, rho ~ p for a perfect gas)
// setExprFields would do the scaling, but v2412 aborts on any expression that reads a field

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "IOmanip.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote("gas mass (sum rho V) of a time, optionally scales p and rho to a given mass");
    argList::addOption("mass", "kg", "gas mass to restore: p and rho are scaled by mass/gasMass");
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

    volScalarField rho
    (
        IOobject("rho", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE),
        mesh
    );
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
