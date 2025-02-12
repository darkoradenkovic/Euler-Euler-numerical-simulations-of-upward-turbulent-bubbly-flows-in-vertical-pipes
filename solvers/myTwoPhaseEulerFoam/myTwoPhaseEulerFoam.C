/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    twoPhaseEulerFoam

Group
    grpMultiphaseSolvers

Description
    Solver for a system of two compressible fluid phases with one dispersed
    phase. Eg, gas bubbles in a liquid.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "myTwoPhaseSystem.H"
#include "PhaseCompressibleMyTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for a system of two incompressible fluid phases with one"
        " dispersed phase.\n"
        "Eg, gas bubbles in a liquid."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"

    bool faceMomentum
    (
        pimple.dict().lookupOrDefault("faceMomentum", false)
    );

    bool implicitPhasePressure
    (
        mesh.solverDict(alpha1.name()).lookupOrDefault
        (
            "implicitPhasePressure", false
        )
    );

    #include "pUf/createDDtU.H"
    #include "pU/createDDtU.H"  
    #include "createDDtUVm.H" 
    #include "calculateVM.H"
   
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve();
            fluid.correct();

            #include "contErrs.H"

            if (faceMomentum)
            {
                #include "pUf/UEqns.H"
                #include "pUf/pEqn.H"
                #include "pUf/DDtU.H"
       
            }
            else
            {
                #include "pU/UEqns.H"
                #include "pU/pEqn.H"
                #include "pU/DDtU.H"
            }

            // calculate material derivatives, to calculate virtual mass force

            #include "DDtUVm.H"
   
         if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
            #include "updateForces.H" 
        }
        
        
        FieldOmega=mag(fvc::grad(phase2.U()));
        FieldS=phase1.d()*FieldOmega/mag(U1-U2);
        URel = U1-U2;
        FiSr=1.0+0.45*FieldS;

        
        
        #include "write.H"

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
