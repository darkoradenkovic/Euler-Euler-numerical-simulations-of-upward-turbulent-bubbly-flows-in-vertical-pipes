/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::mySwarmCorrection

Description

SourceFiles
    mySwarmCorrection.C
    newSwarmCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef mySwarmCorrection_H
#define mySwarmCorrection_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class myPhasePair;

/*---------------------------------------------------------------------------*\
                           Class mySwarmCorrection Declaration
\*---------------------------------------------------------------------------*/

class mySwarmCorrection
{
protected:

    // Protected data

        //- Phase pair
        const myPhasePair& pair_;


public:

    //- Runtime type information
    TypeName("mySwarmCorrection");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        mySwarmCorrection,
        dictionary,
        (
            const dictionary& dict,
            const myPhasePair& pair
        ),
        (dict, pair)
    );


    // Constructors

        //- Construct from a dictionary and a phase pair
        mySwarmCorrection
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    //- Destructor
    virtual ~mySwarmCorrection();


    // Selectors

        static autoPtr<mySwarmCorrection> New
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    // Member Functions

        //- Swarm correction coefficient
        virtual tmp<volScalarField> Cs() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
