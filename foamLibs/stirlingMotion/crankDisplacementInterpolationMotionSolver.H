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

// crank-driven zone interpolation (ESI OpenFOAM v2412)
// bodies move along one axis with a position vs crank angle table (one cycle, deg).
// each time step the gas pressure force on each body (plus the gas drag on porous zones moving with it)
// gives a crank torque F dz/dtheta,
// the flywheel (inertia, load, friction) is stepped and the zones are moved to the new angle.
// points between zones are interpolated linearly along the axis (as displacementInterpolation,
// but zones are matched by name, so the faceZone order in the mesh does not matter)
// remeshing: with stopAngle the run ends with a write when the crank reaches that angle (deg, 0-360),
// the last step lands on it exactly so the next mesh (made at that angle) has the same shape

#ifndef crankDisplacementInterpolationMotionSolver_H
#define crankDisplacementInterpolationMotionSolver_H

#include "displacementMotionSolver.H"
#include "Function1.H"
#include "PtrList.H"
#include "OFstream.H"

namespace Foam
{

class crankDisplacementInterpolationMotionSolver
:
    public displacementMotionSolver
{
    // motion axis (unit)
    vector axis_;

    // bodies: name, patches the gas pushes on, porous zones moving with it (gas drag porDrag_<zone>),
    // displacement vs crank angle
    wordList bodyNames_;
    List<labelList> bodyPatches_;
    List<wordList> bodyPorousZones_;
    PtrList<Function1<scalar>> bodyPosition_;
    scalarField bodyForce_;
    scalarField bodyTorque_;

    // zones: body index (-1 = fixed anchor)
    wordList zoneNames_;
    labelList zoneBody_;

    // point interpolation: breakpoints along the axis (zone -1 = mesh bound)
    labelList rangeZone_;
    labelList pointRange_;
    scalarField pointWeight_;

    // crank state (rad, rad/s)
    scalar theta_;
    scalar omega_;
    scalar thetaOld_;
    scalar omegaOld_;
    label timeIndex_;

    // remeshing: accumulated crank angle to stop at (rad, VGREAT = none) and the time index of the last step
    scalar stopTheta_;
    label stopIndex_;

    // remeshing: the last step of the run is written (restarts are not on the write grid)
    bool endRequested_;
    label endIndex_;

    // flywheel and load
    scalar inertia_;
    word load_;
    scalar omegaFixed_;
    scalar loadTorque_;
    scalar damping_;
    scalar genK_;
    scalar genRLoad_;
    scalar genRInt_;
    scalar frictionTorque_;
    scalar frictionCoeff_;

    // pressure behind the bodies (bounce space, units of p) and density for kinematic p
    scalar pRef_;
    scalar rhoInf_;

    // last step, for the log
    scalar torqueGas_;
    scalar torqueLoad_;
    scalar powerElec_;
    mutable autoPtr<OFstream> logPtr_;

    void readCoeffs();
    void calcInterpolation();
    scalar position(const label bodyi, const scalar theta) const;
    scalar slope(const label bodyi, const scalar theta) const;
    void calcForces();
    void writeLog();
    void writeState() const;

    crankDisplacementInterpolationMotionSolver
    (
        const crankDisplacementInterpolationMotionSolver&
    ) = delete;
    void operator=(const crankDisplacementInterpolationMotionSolver&) = delete;

public:

    TypeName("crankDisplacementInterpolation");

    crankDisplacementInterpolationMotionSolver
    (
        const polyMesh& mesh,
        const IOdictionary& dict
    );

    crankDisplacementInterpolationMotionSolver
    (
        const polyMesh& mesh,
        const IOdictionary& dict,
        const pointVectorField& pointDisplacement,
        const pointIOField& points0
    );

    virtual ~crankDisplacementInterpolationMotionSolver() = default;

    virtual tmp<pointField> curPoints() const;

    virtual void solve();
};

} // End namespace Foam

#endif
