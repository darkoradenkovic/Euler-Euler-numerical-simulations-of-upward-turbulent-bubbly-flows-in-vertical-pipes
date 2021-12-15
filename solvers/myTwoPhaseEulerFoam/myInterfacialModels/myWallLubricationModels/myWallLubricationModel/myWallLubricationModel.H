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
    Foam::myWallLubricationModel

Description

SourceFiles
    myWallLubricationModel.C
    newWallLubricationModel.C

\*---------------------------------------------------------------------------*/

#ifndef myWallLubricationModel_H
#define myWallLubricationModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "myWallDependentModel.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class myPhasePair;

/*---------------------------------------------------------------------------*\
                      Class myWallLubricationModel Declaration
\*---------------------------------------------------------------------------*/

class myWallLubricationModel
:
    public myWallDependentModel
{
protected:

    // Protected data

        //- Phase pair
        const myPhasePair& pair_;


public:

    //- Runtime type information
    TypeName("myWallLubricationModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            myWallLubricationModel,
            dictionary,
            (
                const dictionary& dict,
                const myPhasePair& pair
            ),
            (dict, pair)
        );


    // Static data members

        //- Coefficient dimensions
        static const dimensionSet dimF;


    // Constructors

        //- Construct from components
        myWallLubricationModel
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    //- Destructor
    virtual ~myWallLubricationModel();


    // Selectors

        static autoPtr<myWallLubricationModel> New
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    // Member Functions

        //- Return phase-intensive wall lubrication force
        virtual tmp<volVectorField> Fi() const = 0;

        //- Return wall lubrication force
        virtual tmp<volVectorField> F() const;

        //- Return face wall lubrication force
        virtual tmp<surfaceScalarField> Ff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
