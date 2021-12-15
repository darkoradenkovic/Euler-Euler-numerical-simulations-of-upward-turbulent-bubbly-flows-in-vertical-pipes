volVectorField VirtualMassForce
(
    IOobject
    (
        "VirtualMassForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.Vm()*(DDtU1Vm-DDtU2Vm)
);
