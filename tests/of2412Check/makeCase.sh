#!/bin/bash
# Copyright (c) 2026 Zhi Sheng Chen
#
# Licensed under the PolyForm Noncommercial License 1.0.0 (the "License").
# You may not use this file except in compliance with the License. Any
# commercial use requires a separate commercial license from the copyright
# holder.
#
# You may obtain a copy of the License at:
#     https://polyformproject.org/licenses/noncommercial/1.0.0
#
# Or see the LICENSE file in the root of this repository.
#
# This software is provided "as is", without warranty of any kind, express or
# implied. See the License for the specific language governing permissions
# and limitations.

#writes the dicts for the v2412 check case into the current directory
mkdir -p system constant/tables 0.orig

hdr() { cat <<EOF
FoamFile { version 2.0; format ascii; class $1; object $2; }
EOF
}

#----- controlDict (meshing / motion; solver run overrides application + times)
{ hdr dictionary controlDict; cat <<'EOF'
application     rhoPimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         0.1;
deltaT          0.0025;
writeControl    adjustable;
writeInterval   0.025;
purgeWrite      0;
writeFormat     ascii;
writePrecision  12;
timeFormat      general;
timePrecision   8;
runTimeModifiable false;
adjustTimeStep  no;
maxCo           5;
EOF
} > system/controlDict

{ hdr dictionary blockMeshDict; cat <<'EOF'
vertices ((-0.052 -0.052 -0.002) (0.052 -0.052 -0.002) (0.052 0.052 -0.002) (-0.052 0.052 -0.002)
          (-0.052 -0.052 0.102) (0.052 -0.052 0.102) (0.052 0.052 0.102) (-0.052 0.052 0.102));
blocks (hex (0 1 2 3 4 5 6 7) (26 26 26) simpleGrading (1 1 1));
boundary ( background { type patch; faces ((0 3 2 1) (4 5 6 7) (0 1 5 4) (2 3 7 6) (1 2 6 5) (0 4 7 3)); } );
EOF
} > system/blockMeshDict

{ hdr dictionary decomposeParDict; echo "numberOfSubdomains 4; method scotch;"; } > system/decomposeParDict

#----- snappyHexMeshDict
{ hdr dictionary snappyHexMeshDict; cat <<'EOF'
castellatedMesh true;
snap            true;
addLayers       true;

geometry
{
    WALL-chamber.stl   { type triSurfaceMesh; name WALL-chamber; }
    MOVE-displacer.stl { type triSurfaceMesh; name MOVE-displacer; }
    REGN-matrix.stl    { type triSurfaceMesh; name REGN-matrix; }
    GAPREF-dispSide.stl { type triSurfaceMesh; name GAPREF-dispSide; }
}

castellatedMeshControls
{
    maxLocalCells 2000000;
    maxGlobalCells 20000000;
    minRefinementCells 0;
    maxLoadUnbalance 0.1;
    nCellsBetweenLevels 2;
    features ();
    refinementSurfaces
    {
        WALL-chamber   { level (1 1); patchInfo { type wall; } }
        MOVE-displacer { level (2 2); patchInfo { type wall; } }
        REGN-matrix
        {
            level (2 2);
            faceZone REGN-matrix;     //only used to build the cellZone, not for sampling
            cellZone REGN-matrix;
            cellZoneInside inside;
            faceType internal;
        }
    }
    resolveFeatureAngle 30;
    refinementRegions { GAPREF-dispSide { mode distance; levels ((0.0025 3)); } }
    locationInMesh (0.03 -0.03 0.006);
    allowFreeStandingZoneFaces true;
}

snapControls
{
    nSmoothPatch 3;
    tolerance 2.0;
    nSolveIter 50;
    nRelaxIter 5;
    nFeatureSnapIter 10;
    implicitFeatureSnap true;
    explicitFeatureSnap false;
    multiRegionFeatureSnap false;
}

addLayersControls
{
    relativeSizes true;
    layers
    {
        //precedence test 1: generic first, specific (regex) after
        "WALL-chamber_.*"    { nSurfaceLayers 2; }
        "WALL-chamber_s.*"   { nSurfaceLayers 0; }
        //precedence test 2: specific (regex) first, generic after
        "MOVE-displacer_s.*" { nSurfaceLayers 0; }
        "MOVE-displacer_.*"  { nSurfaceLayers 2; }
    }
    expansionRatio 1.2;
    finalLayerThickness 0.4;
    minThickness 0.1;
    nGrow 0;
    featureAngle 60;
    slipFeatureAngle 30;
    nRelaxIter 3;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedialAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
}

meshQualityControls
{
    #includeEtc "caseDicts/meshQualityDict"
    nSmoothScale 4;
    errorReduction 0.75;
}

mergeTolerance 1e-6;
EOF
} > system/snappyHexMeshDict

#----- topoSetDict: motion zones from patches (regex), plane + surface sample zones
{ hdr dictionary topoSetDict; cat <<'EOF'
actions
(
    { name hotSet;   type faceSet; action new; source patchToFace; patch "WALL-chamber_h.*"; }
    { name hotAnchor; type faceZoneSet; action new; source setToFaceZone; faceSet hotSet; }
    { name coldSet;  type faceSet; action new; source patchToFace; patch "WALL-chamber_c[o].*"; }
    { name coldAnchor; type faceZoneSet; action new; source setToFaceZone; faceSet coldSet; }
    { name crownSet; type faceSet; action new; source patchToFace; patch "WALL-chamber_crown"; }
    { name crown; type faceZoneSet; action new; source setToFaceZone; faceSet crownSet; }
    { name dTopSet;  type faceSet; action new; source patchToFace; patch "MOVE-displacer_t.*"; }
    { name displacerTop; type faceZoneSet; action new; source setToFaceZone; faceSet dTopSet; }
    { name dBotSet;  type faceSet; action new; source patchToFace; patch "MOVE-displacer_b.*"; }
    { name displacerBottom; type faceZoneSet; action new; source setToFaceZone; faceSet dBotSet; }

    { name planeMid; type faceZoneSet; action new; source planeToFaceZone;
      point (0 0 0.0301); normal (0 0 1); option all; }
    { name smpMid; type faceZoneSet; action new; source searchableSurfaceToFaceZone;
      surfaceType triSurfaceMesh; surfaceName SMP-mid.stl; }

    //regenerator outlet sample zones from the per-PID cap STLs (select mode)
    { name capCold; type faceZoneSet; action new; source searchableSurfaceToFaceZone;
      surfaceType triSurfaceMesh; surfaceName SMP-capCold.stl; }
    { name capHot; type faceZoneSet; action new; source searchableSurfaceToFaceZone;
      surfaceType triSurfaceMesh; surfaceName SMP-capHot.stl; }
);
EOF
} > system/topoSetDict

#----- mesh motion
{ hdr dictionary dynamicMeshDict; cat <<'EOF'
dynamicFvMesh   dynamicMotionSolverFvMesh;
motionSolver    displacementInterpolation;
displacementInterpolationCoeffs
{
    interpolationTables
    (
        (hotAnchor       zeroTable)
        (coldAnchor      zeroTable)
        (displacerTop    displacerTable)
        (displacerBottom displacerTable)
        (crown           pistonTable)
    );
}
EOF
} > constant/dynamicMeshDict

#----- thermo / transport
{ hdr dictionary thermophysicalProperties; cat <<'EOF'
thermoType { type hePsiThermo; mixture pureMixture; transport const; thermo hConst;
             equationOfState perfectGas; specie specie; energy sensibleInternalEnergy; }
mixture
{
    specie         { molWeight 28.96; }
    thermodynamics { Cp 1005; Hf 0; }
    transport      { mu 1.8e-05; Pr 0.7; }
}
EOF
} > constant/thermophysicalProperties
{ hdr dictionary turbulenceProperties; echo "simulationType laminar;"; } > constant/turbulenceProperties

#----- fvOptions: porosity + LTNE coded source in the regenerator zone
{ hdr dictionary fvOptions; cat <<'EOF'
regenPorosity
{
    type explicitPorositySource;
    explicitPorositySourceCoeffs
    {
        selectionMode cellZone;
        cellZone REGN-matrix;
        type DarcyForchheimer;
        DarcyForchheimerCoeffs
        {
            d (1e8 1e8 1e8);
            f (1e4 1e4 1e4);
            coordinateSystem { origin (0 0 0); e1 (1 0 0); e2 (0 1 0); }
        }
    }
}

regenLTNE
{
    type scalarCodedSource;
    name regenLTNE;
    scalarCodedSourceCoeffs
    {
        selectionMode cellZone;
        cellZone REGN-matrix;
        fields (e);

        codeInclude
        #{
            #include "fvCFD.H"
        #};

        codeCorrect #{ #};
        codeAddSup #{ #};
        codeConstrain #{ #};

        codeAddSupRho
        #{
            //matrix temperature Ts kept as a registered field between steps
            const fvMesh& m = mesh_;
            const volScalarField& T = m.lookupObject<volScalarField>("T");
            if (!m.foundObject<volScalarField>("Ts"))
            {
                volScalarField* TsPtr = new volScalarField
                (
                    IOobject("Ts", m.time().timeName(), m,
                        IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                    T
                );
                TsPtr->store();
            }
            volScalarField& Ts = const_cast<volScalarField&>(m.lookupObject<volScalarField>("Ts"));

            const scalar ha = 1e6;          //h*a [W/m3/K]
            const scalar Cs = 0.1*7900*500; //(1-eps)*rhos*cs [J/m3/K]
            const scalar dt = m.time().deltaTValue();
            const labelList& cells = this->cells();
            const scalarField& V = m.V();
            scalarField& src = eqn.source();
            forAll(cells, i)
            {
                const label c = cells[i];
                const scalar Tg = T[c];
                //exact update of dTs/dt = ha/Cs*(Tg - Ts)
                Ts[c] = Tg + (Ts[c] - Tg)*Foam::exp(-ha*dt/Cs);
                src[c] -= ha*(Ts[c] - Tg)*V[c];
            }
            Info<< "regenLTNE: mean Ts " << gAverage(scalarField(Ts, cells)) << endl;
        #};
    }
}
EOF
} > constant/fvOptions

#----- numerics
{ hdr dictionary fvSchemes; cat <<'EOF'
ddtSchemes      { default Euler; }
gradSchemes     { default Gauss linear; }
divSchemes
{
    default none;
    div(phi,U) Gauss limitedLinearV 1;
    div(phi,e) Gauss limitedLinear 1;
    div(phi,K) Gauss linear;
    div(phiv,p) Gauss linear;
    div(meshPhi,p) Gauss linear;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes   { default limited 0.5; }
EOF
} > system/fvSchemes
{ hdr dictionary fvSolution; cat <<'EOF'
solvers
{
    "rho.*" { solver diagonal; }
    p      { solver GAMG; smoother GaussSeidel; tolerance 1e-2; relTol 0.01; }
    pFinal { $p; relTol 0; }
    "(U|e).*" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-8; relTol 0; }
}
PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 1;
    nCorrectors 2;
    nNonOrthogonalCorrectors 1;
    correctPhi no;
}
EOF
} > system/fvSolution

#----- initial fields
{ hdr volScalarField p; cat <<'EOF'
dimensions [1 -1 -2 0 0 0 0];
internalField uniform 1e5;
boundaryField { ".*" { type zeroGradient; } }
EOF
} > 0.orig/p
{ hdr volScalarField T; cat <<'EOF'
dimensions [0 0 0 1 0 0 0];
internalField uniform 300;
boundaryField
{
    ".*" { type zeroGradient; }
    "WALL-chamber_hot"  { type fixedValue; value uniform 600; }
    "WALL-chamber_cold" { type fixedValue; value uniform 300; }
    "WALL-chamber_liner"
    {
        type externalWallHeatFluxTemperature;
        mode coefficient;
        Ta constant 300;
        h constant 15;
        kappaMethod fluidThermo;
        value $internalField;
    }
}
EOF
} > 0.orig/T
{ hdr volVectorField U; cat <<'EOF'
dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
    ".*" { type noSlip; }
    "WALL-chamber_crown" { type movingWallVelocity; value uniform (0 0 0); }
    "MOVE-displacer_.*"  { type movingWallVelocity; value uniform (0 0 0); }
}
EOF
} > 0.orig/U

#----- run-time sampling on faceZones (select-mode zones only)
cat > system/samplingFOs <<'EOF'
capColdFlow  { type surfaceFieldValue; libs (fieldFunctionObjects); regionType faceZone; name capCold;
               operation sum; fields (phi); writeFields false; writeControl timeStep; writeInterval 1; log true; }
capHotFlow   { $capColdFlow; name capHot; }
smpMidFlow   { $capColdFlow; name smpMid; }
planeMidFlow { $capColdFlow; name planeMid; }
capColdAvg   { type surfaceFieldValue; libs (fieldFunctionObjects); regionType faceZone; name capCold;
               operation areaAverage; fields (rho T p); writeFields false; writeControl timeStep; writeInterval 1; log true; }
gasMass      { type volFieldValue; libs (fieldFunctionObjects); regionType all;
               operation volIntegrate; fields (rho); writeFields false; writeControl timeStep; writeInterval 1; log true; }
EOF
#----- meshPhi sums: meshPhi is not visible to function objects at run time (v2412),
#      so these run with postProcess -fields '(phi meshPhi)' on written times
cat > system/meshFluxFOs <<'EOF'
FoamFile { version 2.0; format ascii; class dictionary; object meshFluxFOs; }
functions
{
    capColdMesh  { type surfaceFieldValue; libs (fieldFunctionObjects); regionType faceZone; name capCold;
                   operation sum; fields (meshPhi phi); writeFields false; log true; }
    capHotMesh   { $capColdMesh; name capHot; }
    smpMidMesh   { $capColdMesh; name smpMid; }
    planeMidMesh { $capColdMesh; name planeMid; }
}
EOF
echo "case dicts written"
