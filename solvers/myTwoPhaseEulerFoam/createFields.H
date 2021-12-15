#include "readGravitationalAcceleration.H"
#include "readhRef.H"

Info<< "Creating myTwoPhaseSystem\n" << endl;

myTwoPhaseSystem fluid(mesh, g);

myPhaseModel& phase1 = fluid.phase1();
myPhaseModel& phase2 = fluid.phase2();

volScalarField& alpha1 = phase1;
volVectorField& U1 = phase1.U();
volVectorField& U2 = phase2.U();

volScalarField& p = phase1.thermo().p();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    fluid
);

#include "gh.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



volScalarField FieldOmega
(
    IOobject
    (
        "FieldOmega",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(fvc::grad(U2))
);

volScalarField dB
(
    IOobject
    (
        "dB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phase1.d()
);


volScalarField FieldS
(
    IOobject
    (
        "FieldS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    phase1.d()*FieldOmega/mag(U1-U2)
);

volVectorField URel
(
    IOobject
    (
        "URel",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U1-U2
);

volScalarField FiSr
(
    IOobject
    (
        "FiSr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
   1.0+0.45*FieldS
);


volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.U()
);


volVectorField DragForce
(
    IOobject
    (
        "DragForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.Kd()*(U1-U2)
);

volVectorField LiftForce
(
    IOobject
    (
        "LiftForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.LiftForce()
);

volVectorField WallLubricationForce
(
    IOobject
    (
        "WallLubricationForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.WallLubricationForce()
);

volVectorField LamPlusTurbDispForce
(
    IOobject
    (
        "LamPlusTurbDispForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.LamPlusTurbDispForce()
);


volVectorField LamDispForce
(
    IOobject
    (
        "LamDispForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.LamDispForce()
);

volVectorField TurbDispForce
(
    IOobject
    (
        "TurbDispForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.TurbDispForce()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);
mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(p.dimensions()/dimTime, Zero)
);


Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K1(IOobject::groupName("K", phase1.name()), 0.5*magSqr(U1));
volScalarField K2(IOobject::groupName("K", phase2.name()), 0.5*magSqr(U2));

#include "createMRF.H"
#include "createFvOptions.H"
