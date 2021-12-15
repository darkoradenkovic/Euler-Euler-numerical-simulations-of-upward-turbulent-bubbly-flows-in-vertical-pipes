
    volVectorField DDtU1Vm
    (
        "DDtU1Vm",
        fvc::ddt(U1)
      + fvc::div(phi1, U1)
      - fvc::div(phi1)*U1
    );

    volVectorField DDtU2Vm
    (
        "DDtU2Vm",
        fvc::ddt(U2)
      + fvc::div(phi2, U2)
      - fvc::div(phi2)*U2
    );
