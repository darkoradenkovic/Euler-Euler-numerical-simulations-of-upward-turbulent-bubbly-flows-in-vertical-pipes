/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Class
    Foam::myDiameterModel

Description
    A2stract base-class for dispersed-phase particle diameter models.

SourceFiles
    myDiameterModel.C
    newDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef myDiameterModel_H
#define myDiameterModel_H

#include "dictionary.H"
#include "myPhaseModel.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class myDiameterModel Declaration
\*---------------------------------------------------------------------------*/

class myDiameterModel
{

protected:

    // Protected data

        dictionary diameterProperties_;
        const myPhaseModel& phase_;


public:

    //- Runtime type information
    TypeName("myDiameterModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            myDiameterModel,
            dictionary,
            (
                const dictionary& diameterProperties,
                const myPhaseModel& phase
            ),
            (diameterProperties, phase)
        );


    // Constructors

        myDiameterModel
        (
            const dictionary& diameterProperties,
            const myPhaseModel& phase
        );


    //- Destructor
    virtual ~myDiameterModel();


    // Selectors

        static autoPtr<myDiameterModel> New
        (
            const dictionary& diameterProperties,
            const myPhaseModel& phase
        );


    // Member Functions

        //- Return the phase diameter properties dictionary
        const dictionary& diameterProperties() const
        {
            return diameterProperties_;
        }

        //- Return the phase
        const myPhaseModel& phase() const
        {
            return phase_;
        }

        //- Return the phase mean diameter field
        virtual tmp<volScalarField> d() const = 0;

        //- Correct the diameter field
        virtual void correct();

        //- Read phaseProperties dictionary
        virtual bool read(const dictionary& phaseProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
