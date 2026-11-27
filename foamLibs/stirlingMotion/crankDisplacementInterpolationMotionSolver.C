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

#include "crankDisplacementInterpolationMotionSolver.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"
#include "volFields.H"
#include "unitConversion.H"
#include "OSspecific.H"

namespace Foam
{
    defineTypeNameAndDebug(crankDisplacementInterpolationMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        crankDisplacementInterpolationMotionSolver,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        displacementMotionSolver,
        crankDisplacementInterpolationMotionSolver,
        displacement
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::crankDisplacementInterpolationMotionSolver::readCoeffs()
{
    const dictionary& coeffs = coeffDict();

    axis_ = coeffs.getOrDefault<vector>("axis", vector(0, 0, 1));
    axis_.normalise();

    // bodies
    const dictionary& bodiesDict = coeffs.subDict("bodies");
    bodyNames_ = bodiesDict.toc();
    bodyPatches_.setSize(bodyNames_.size());
    bodyPorousZones_.setSize(bodyNames_.size());
    bodyPosition_.setSize(bodyNames_.size());
    bodyForce_.setSize(bodyNames_.size(), Zero);
    bodyTorque_.setSize(bodyNames_.size(), Zero);

    forAll(bodyNames_, bodyi)
    {
        const dictionary& bodyDict = bodiesDict.subDict(bodyNames_[bodyi]);
        const wordRes patchNames(bodyDict.get<wordRes>("patches"));
        bodyPatches_[bodyi] = mesh().boundaryMesh().indices(patchNames);
        if (bodyPatches_[bodyi].empty())
        {
            FatalIOErrorInFunction(bodyDict)
                << "body " << bodyNames_[bodyi] << ": no patches match "
                << patchNames << nl << "Known patches: "
                << mesh().boundaryMesh().names() << exit(FatalIOError);
        }
        bodyPorousZones_[bodyi] = bodyDict.getOrDefault<wordList>("porousZones", wordList());
        for (const word& zoneName : bodyPorousZones_[bodyi])
        {
            if (mesh().cellZones().findZoneID(zoneName) == -1)
            {
                FatalIOErrorInFunction(bodyDict)
                    << "body " << bodyNames_[bodyi] << ": no cellZone " << zoneName << nl
                    << "Known cellZones: " << mesh().cellZones().names() << exit(FatalIOError);
            }
        }
        bodyPosition_.set
        (
            bodyi,
            Function1<scalar>::New("position", bodyDict, &mesh())
        );
    }

    // zones
    const List<Pair<word>> zoneList(coeffs.lookup("zones"));
    zoneNames_.setSize(zoneList.size());
    zoneBody_.setSize(zoneList.size());
    forAll(zoneList, zonei)
    {
        zoneNames_[zonei] = zoneList[zonei].first();
        const word& bodyName = zoneList[zonei].second();
        zoneBody_[zonei] = (bodyName == "none" ? -1 : bodyNames_.find(bodyName));
        if (bodyName != "none" && zoneBody_[zonei] == -1)
        {
            FatalIOErrorInFunction(coeffs)
                << "zone " << zoneNames_[zonei] << ": body " << bodyName
                << " not in bodies " << bodyNames_ << exit(FatalIOError);
        }
    }

    // flywheel and load
    const dictionary& crankDict = coeffs.subDict("crank");
    theta_ = degToRad(crankDict.get<scalar>("theta0"));
    omega_ = crankDict.get<scalar>("omega0");
    inertia_ = crankDict.get<scalar>("inertia");
    load_ = crankDict.get<word>("load");
    omegaFixed_ = omega_;
    loadTorque_ = crankDict.getOrDefault<scalar>("torque", 0);
    damping_ = crankDict.getOrDefault<scalar>("damping", 0);
    frictionTorque_ = crankDict.getOrDefault<scalar>("frictionTorque", 0);
    frictionCoeff_ = crankDict.getOrDefault<scalar>("frictionCoeff", 0);
    pRef_ = crankDict.getOrDefault<scalar>("pRef", 0);
    rhoInf_ = crankDict.getOrDefault<scalar>("rhoInf", 1);
    genK_ = genRLoad_ = genRInt_ = 0;

    if (load_ == "generator")
    {
        const dictionary& genDict = crankDict.subDict("generator");
        genK_ = genDict.get<scalar>("k");
        genRLoad_ = genDict.get<scalar>("RLoad");
        genRInt_ = genDict.getOrDefault<scalar>("RInternal", 0);
        if (genRLoad_ + genRInt_ <= 0)
        {
            FatalIOErrorInFunction(genDict)
                << "RLoad + RInternal must be > 0" << exit(FatalIOError);
        }
    }
    else if (load_ != "speed" && load_ != "torque" && load_ != "damper")
    {
        FatalIOErrorInFunction(crankDict)
            << "load " << load_ << " unknown, use speed, torque, damper or generator"
            << exit(FatalIOError);
    }
    if (load_ != "speed" && inertia_ <= 0)
    {
        FatalIOErrorInFunction(crankDict)
            << "inertia must be > 0" << exit(FatalIOError);
    }

    // restart: crank state of the start time
    IOdictionary state
    (
        IOobject
        (
            "crankState",
            mesh().time().timeName(),
            "uniform",
            mesh().time(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false    // not registered (bool: IOobject::NO_REGISTER only exists from v2306)
        )
    );
    if (state.found("theta"))
    {
        theta_ = degToRad(state.get<scalar>("theta"));
        omega_ = state.get<scalar>("omega");
        Info<< "crankDisplacementInterpolation: restart at crank angle "
            << radToDeg(theta_) << " deg, omega " << omega_ << " rad/s" << endl;
    }
    thetaOld_ = theta_;
    omegaOld_ = omega_;
    timeIndex_ = mesh().time().timeIndex();
    torqueGas_ = torqueLoad_ = powerElec_ = 0;

    // remeshing: the next time the crank reaches stopAngle, at least a millidegree ahead
    stopTheta_ = VGREAT;
    stopIndex_ = labelMax;
    endRequested_ = false;
    endIndex_ = labelMax;
    if (crankDict.found("stopAngle"))
    {
        const scalar stopDeg = crankDict.get<scalar>("stopAngle");
        const scalar thetaDeg = radToDeg(theta_);
        stopTheta_ = degToRad(stopDeg + 360*std::ceil((thetaDeg + 1e-3 - stopDeg)/360));
        Info<< "crankDisplacementInterpolation: stops at crank angle " << radToDeg(stopTheta_)
            << " deg (remeshing)" << endl;
    }
}


void Foam::crankDisplacementInterpolationMotionSolver::calcInterpolation()
{
    const faceZoneMesh& fZones = mesh().faceZones();
    const scalarField coord(points0() & axis_);

    // two breakpoints per zone (its min and max along the axis)
    SortableList<scalar> zoneCoord(2*zoneNames_.size());
    forAll(zoneNames_, zonei)
    {
        const label zoneID = fZones.findZoneID(zoneNames_[zonei]);
        if (zoneID == -1)
        {
            FatalErrorInFunction
                << "Cannot find zone " << zoneNames_[zonei] << nl
                << "Valid zones are " << fZones.names() << exit(FatalError);
        }
        scalar minCoord = VGREAT;
        scalar maxCoord = -VGREAT;
        for (const label pointi : fZones[zoneID]().meshPoints())
        {
            minCoord = min(minCoord, coord[pointi]);
            maxCoord = max(maxCoord, coord[pointi]);
        }
        zoneCoord[2*zonei] = returnReduce(minCoord, minOp<scalar>());
        zoneCoord[2*zonei + 1] = returnReduce(maxCoord, maxOp<scalar>());
    }
    zoneCoord.sort();
    zoneCoord.first() -= SMALL;
    zoneCoord.last() += SMALL;

    // add the mesh bounds (zero displacement) outside the outermost zones
    const scalar minMesh = gMin(coord) - SMALL;
    const scalar maxMesh = gMax(coord) + SMALL;
    DynamicList<scalar> rangeCoord(zoneCoord.size() + 2);
    DynamicList<label> rangeZone(zoneCoord.size() + 2);
    if (minMesh < zoneCoord.first())
    {
        rangeCoord.append(minMesh);
        rangeZone.append(-1);
    }
    forAll(zoneCoord, i)
    {
        rangeCoord.append(zoneCoord[i]);
        rangeZone.append(zoneCoord.indices()[i]/2);
    }
    if (maxMesh > zoneCoord.last())
    {
        rangeCoord.append(maxMesh);
        rangeZone.append(-1);
    }
    rangeZone_.transfer(rangeZone);

    pointRange_.setSize(coord.size());
    pointWeight_.setSize(coord.size());
    forAll(coord, pointi)
    {
        const label rangei = findLower(rangeCoord, coord[pointi]);
        if (rangei == -1 || rangei == rangeCoord.size() - 1)
        {
            FatalErrorInFunction
                << "point " << points0()[pointi] << " outside the zone ranges "
                << rangeCoord << abort(FatalError);
        }
        pointRange_[pointi] = rangei;
        pointWeight_[pointi] =
            (coord[pointi] - rangeCoord[rangei])
           /(rangeCoord[rangei + 1] - rangeCoord[rangei]);
    }

    Info<< "crankDisplacementInterpolation: " << zoneNames_.size() << " zones, "
        << bodyNames_.size() << " bodies, load " << load_ << endl;
}


Foam::scalar Foam::crankDisplacementInterpolationMotionSolver::position
(
    const label bodyi,
    const scalar theta
) const
{
    // tables are in deg over one cycle (outOfBounds repeat)
    return bodyPosition_[bodyi].value(radToDeg(theta));
}


Foam::scalar Foam::crankDisplacementInterpolationMotionSolver::slope
(
    const label bodyi,
    const scalar theta
) const
{
    // dz/dtheta (m/rad), central difference over half a degree
    const scalar h = degToRad(0.25);
    return (position(bodyi, theta + h) - position(bodyi, theta - h))/(2*h);
}


void Foam::crankDisplacementInterpolationMotionSolver::calcForces()
{
    // pressure force of the gas on each body along the axis (boundary Sf points into the body)
    bodyForce_ = Zero;
    const volScalarField* pPtr = mesh().findObject<volScalarField>("p");
    if (!pPtr)
    {
        return;
    }
    const volScalarField& p = *pPtr;
    const scalar factor = (p.dimensions() == dimPressure ? 1.0 : rhoInf_);

    forAll(bodyPatches_, bodyi)
    {
        scalar force = 0;
        for (const label patchi : bodyPatches_[bodyi])
        {
            const vectorField Sf(mesh().boundaryMesh()[patchi].faceAreas());
            force += sum((p.boundaryField()[patchi] - pRef_)*(Sf & axis_));
        }
        // gas drag on porous zones moving with the body (N/m3, from the zone's drag source, last solved U)
        for (const word& zoneName : bodyPorousZones_[bodyi])
        {
            const volVectorField* dragPtr = mesh().findObject<volVectorField>("porDrag_" + zoneName);
            if (!dragPtr)
            {
                continue;
            }
            const scalarField& V = refCast<const fvMesh>(mesh()).V();
            for (const label celli : mesh().cellZones()[mesh().cellZones().findZoneID(zoneName)])
            {
                force += ((*dragPtr)[celli] & axis_)*V[celli]/factor;
            }
        }
        bodyForce_[bodyi] = factor*returnReduce(force, sumOp<scalar>());
    }
}


void Foam::crankDisplacementInterpolationMotionSolver::writeLog()
{
    if (!Pstream::master())
    {
        return;
    }
    if (!logPtr_)
    {
        const fileName dir
        (
            mesh().time().globalPath()/"postProcessing"/"crank"
           /mesh().time().timeName(mesh().time().startTime().value())
        );
        mkDir(dir);
        logPtr_.reset(new OFstream(dir/"crank.dat"));
        OFstream& os = *logPtr_;
        os  << "# crank angle (deg, accumulated), omega (rad/s), torques (N m, gas + = drives the crank),"
            << " electrical power (W), body forces along the axis (N), body torques on the crank (N m)" << nl
            << "# Time\ttheta\tomega\ttorqueGas\ttorqueLoad\tpowerElec";
        forAll(bodyNames_, bodyi)
        {
            os  << "\tF_" << bodyNames_[bodyi];
        }
        forAll(bodyNames_, bodyi)
        {
            os  << "\tT_" << bodyNames_[bodyi];
        }
        os  << nl;
    }
    OFstream& os = *logPtr_;
    os  << mesh().time().value() << tab << radToDeg(theta_) << tab << omega_
        << tab << torqueGas_ << tab << torqueLoad_ << tab << powerElec_;
    forAll(bodyForce_, bodyi)
    {
        os  << tab << bodyForce_[bodyi];
    }
    forAll(bodyTorque_, bodyi)
    {
        os  << tab << bodyTorque_[bodyi];
    }
    os  << endl;
}


void Foam::crankDisplacementInterpolationMotionSolver::writeState() const
{
    IOdictionary state
    (
        IOobject
        (
            "crankState",
            mesh().time().timeName(),
            "uniform",
            mesh().time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false    // not registered
        )
    );
    state.add("theta", radToDeg(theta_));
    state.add("omega", omega_);
    state.regIOobject::write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::crankDisplacementInterpolationMotionSolver::
crankDisplacementInterpolationMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName)
{
    readCoeffs();
    calcInterpolation();
}


Foam::crankDisplacementInterpolationMotionSolver::
crankDisplacementInterpolationMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict,
    const pointVectorField& pointDisplacement,
    const pointIOField& points0
)
:
    displacementMotionSolver(mesh, dict, pointDisplacement, points0, typeName)
{
    readCoeffs();
    calcInterpolation();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::crankDisplacementInterpolationMotionSolver::curPoints() const
{
    if (mesh().nPoints() != points0().size())
    {
        FatalErrorInFunction
            << "The number of points in the mesh has changed." << exit(FatalError);
    }

    scalarField zoneDisp(zoneBody_.size(), Zero);
    forAll(zoneBody_, zonei)
    {
        if (zoneBody_[zonei] >= 0)
        {
            zoneDisp[zonei] = position(zoneBody_[zonei], theta_);
        }
    }

    auto tcurPoints = tmp<pointField>::New(points0());
    auto& curPoints = tcurPoints.ref();
    forAll(curPoints, pointi)
    {
        const label rangei = pointRange_[pointi];
        const label zoneLo = rangeZone_[rangei];
        const label zoneHi = rangeZone_[rangei + 1];
        const scalar dLo = (zoneLo == -1 ? 0.0 : zoneDisp[zoneLo]);
        const scalar dHi = (zoneHi == -1 ? 0.0 : zoneDisp[zoneHi]);
        const scalar w = pointWeight_[pointi];
        curPoints[pointi] += ((1.0 - w)*dLo + w*dHi)*axis_;
    }
    return tcurPoints;
}


void Foam::crankDisplacementInterpolationMotionSolver::solve()
{
    // called at the start of each time step (again with moveMeshOuterCorrectors:
    // then the step is redone from the old state with the newer pressure)
    const Time& runTime = mesh().time();
    if (runTime.timeIndex() != timeIndex_)
    {
        thetaOld_ = theta_;
        omegaOld_ = omega_;
        timeIndex_ = runTime.timeIndex();
    }
    const scalar dt = runTime.deltaTValue();

    // gas torque on the crank from the current pressure (mesh at the old angle)
    calcForces();
    forAll(bodyForce_, bodyi)
    {
        bodyTorque_[bodyi] = bodyForce_[bodyi]*slope(bodyi, thetaOld_);
    }
    torqueGas_ = sum(bodyTorque_);
    const scalar friction = frictionTorque_*sign(omegaOld_);

    // flywheel over a step dtStep: I domega/dt = torqueGas - load - friction, speed-proportional parts implicit
    auto advance = [&](const scalar dtStep)
    {
        powerElec_ = 0;
        if (load_ == "speed")
        {
            omega_ = omegaFixed_;
            torqueLoad_ = torqueGas_;
        }
        else
        {
            scalar c = frictionCoeff_;
            scalar T0 = friction;
            if (load_ == "torque")
            {
                T0 += loadTorque_;
            }
            else if (load_ == "damper")
            {
                c += damping_;
            }
            else
            {
                c += sqr(genK_)/(genRLoad_ + genRInt_);
            }
            omega_ = (inertia_*omegaOld_ + dtStep*(torqueGas_ - T0))/(inertia_ + dtStep*c);
            torqueLoad_ = T0 + c*omega_;
            if (load_ == "generator")
            {
                const scalar current = genK_*omega_/(genRLoad_ + genRInt_);
                powerElec_ = sqr(current)*genRLoad_;
            }
        }
        theta_ = thetaOld_ + 0.5*dtStep*(omegaOld_ + omega_);
    };
    advance(dt);

    // the last step of a segment or of the run is shortened (or lengthened, up to 1.2x) to end exactly on the stop
    // angle or endTime, so time and crank angle stay together (with a fixed speed the angle is omega t: setting only
    // the angle lost 0.2 deg per map). Called before the equations of the step are solved
    auto land = [&](const scalar dtNew)
    {
        if (dtNew > SMALL && mag(dtNew - dt) > 1e-12*dt)
        {
            Time& time = const_cast<Time&>(runTime);
            time.setTime(runTime.value() - dt + dtNew, runTime.timeIndex());
            time.setDeltaT(dtNew, false);
            advance(dtNew);
        }
    };

    // remeshing: the step after the stop is requested is the last one and lands exactly on the stop angle
    // (it moves up to 1.5 normal steps), the stop is requested when the next step could reach it
    if (stopTheta_ < 0.5*VGREAT)
    {
        if (runTime.timeIndex() >= stopIndex_)
        {
            //step length that reaches the stop angle (exact for a fixed speed, close for the flywheel)
            land(dt*(stopTheta_ - thetaOld_)/max(theta_ - thetaOld_, SMALL));
            theta_ = stopTheta_;
        }
        else if (theta_ >= stopTheta_)
        {
            //passed before the stop was requested (the time step stretched to a write time can grow more than
            //1.5x): go half way now and land on the stop angle next step, the stop only takes effect after the
            //next step, so landing now gave an extra step with the crank at rest
            land(0.5*dt*(stopTheta_ - thetaOld_)/max(theta_ - thetaOld_, SMALL));
        }
        else if (runTime.timeIndex() >= endIndex_)
        {
            land(runTime.endTime().value() - (runTime.value() - dt));
        }
        //1.5 steps ahead: the time step grows by at most 1.2x per step, so the stop angle is normally not passed
        //before the stop is requested (write-time stretching can still do it, handled above)
        if (stopIndex_ == labelMax && theta_ + 1.5*(theta_ - thetaOld_) >= stopTheta_)
        {
            stopIndex_ = runTime.timeIndex() + 1;
            runTime.stopAt(Time::stopAtControls::saWriteNow);
            Info<< "crankDisplacementInterpolation: next step reaches the stop angle "
                << radToDeg(stopTheta_) << " deg, writing and stopping after it" << endl;
        }

        // end of the run: written too, the write grid starts at each restart time so endTime is rarely on it
        const scalar dtNow = runTime.deltaTValue();
        if
        (
            !endRequested_ && stopIndex_ == labelMax
         && runTime.value() + 1.7*dtNow >= runTime.endTime().value()
        )
        {
            endRequested_ = true;
            endIndex_ = runTime.timeIndex() + 1;
            runTime.stopAt(Time::stopAtControls::saWriteNow);
            Info<< "crankDisplacementInterpolation: next step reaches endTime, writing and stopping after it"
                << endl;
        }
    }

    writeLog();
    if (runTime.writeTime())
    {
        writeState();
    }
}


// ************************************************************************* //
