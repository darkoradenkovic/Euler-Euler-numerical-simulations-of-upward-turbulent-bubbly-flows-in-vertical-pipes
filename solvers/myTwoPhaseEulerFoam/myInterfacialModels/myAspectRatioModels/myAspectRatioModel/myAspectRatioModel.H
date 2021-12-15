/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2015 OpenFOAM Foundation
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
    Foam::myAspectRatioModel

Description

SourceFiles
    myAspectRatioModel.C
    newAspectRatioModel.C

\*---------------------------------------------------------------------------*/

#ifndef myAspectRatioModel_H
#define myAspectRatioModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class myPhasePair;

/*---------------------------------------------------------------------------*\
                           Class myAspectRatioModel Declaration
\*---------------------------------------------------------------------------*/

class myAspectRatioModel
{
protected:

    // Protected data

        //- Phase pair
        const myPhasePair& pair_;


public:

    //- Runtime type information
    TypeName("myAspectRatioModel");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        myAspectRatioModel,
        dictionary,
        (
            const dictionary& dict,
            const myPhasePair& pair
        ),
        (dict, pair)
    );


    // Constructors

        //- Construct from a dictionary and an ordered phase pair
        myAspectRatioModel
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    //- Destructor
    virtual ~myAspectRatioModel();


    // Selectors

        static autoPtr<myAspectRatioModel> New
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    // Member Functions

        //- Aspect ratio
        virtual tmp<volScalarField> E() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
