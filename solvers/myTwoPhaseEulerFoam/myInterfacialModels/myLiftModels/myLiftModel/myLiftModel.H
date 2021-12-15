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
    Foam::myLiftModel

Description

SourceFiles
    myLiftModel.C
    newLiftModel.C

\*---------------------------------------------------------------------------*/

#ifndef myLiftModel_H
#define myLiftModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class myPhasePair;

/*---------------------------------------------------------------------------*\
                           Class myLiftModel Declaration
\*---------------------------------------------------------------------------*/

class myLiftModel
{
protected:

    // Protected data

        //- Phase pair
        const myPhasePair& pair_;


public:

    //- Runtime type information
    TypeName("myLiftModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            myLiftModel,
            dictionary,
            (
                const dictionary& dict,
                const myPhasePair& pair
            ),
            (dict, pair)
        );


    // Static data members

        //- Force dimensions
        static const dimensionSet dimF;


    // Constructors

        //- Construct from a dictionary and a phase pair
        myLiftModel
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    //- Destructor
    virtual ~myLiftModel();


    // Selectors

        static autoPtr<myLiftModel> New
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    // Member Functions

        //- Return lift coefficient
        virtual tmp<volScalarField> Cl() const = 0;

        //- Return phase-intensive lift force
        virtual tmp<volVectorField> Fi() const;

        //- Return lift force
        virtual tmp<volVectorField> F() const;

        //- Return face lift force
        virtual tmp<surfaceScalarField> Ff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
