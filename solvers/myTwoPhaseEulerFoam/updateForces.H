    DragForce=fluid.Kd()*(U1-U2);
    LiftForce=fluid.LiftForce();
    WallLubricationForce=fluid.WallLubricationForce();
    LamPlusTurbDispForce=fluid.LamPlusTurbDispForce();
    LamDispForce=fluid.LamDispForce();
    TurbDispForce=fluid.TurbDispForce();
    VirtualMassForce=fluid.Vm()*(DDtU1Vm-DDtU2Vm);
