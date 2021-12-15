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
    Foam::myTurbulentDispersionModel

Description

SourceFiles
    myTurbulentDispersionModel.C
    newTurbulentDispersionModel.C

\*---------------------------------------------------------------------------*/

#ifndef myTurbulentDispersionModel_H
#define myTurbulentDispersionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class myPhasePair;

/*---------------------------------------------------------------------------*\
                    Class myTurbulentDispersionModel Declaration
\*---------------------------------------------------------------------------*/

class myTurbulentDispersionModel
{
protected:

    // Protected data

        //- Phase pair
        const myPhasePair& pair_;


public:

    //- Runtime type information
    TypeName("myTurbulentDispersionModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            myTurbulentDispersionModel,
            dictionary,
            (
                const dictionary& dict,
                const myPhasePair& pair
            ),
            (dict, pair)
        );

    // Static data members

        //- Diffusivity dimensions
        static const dimensionSet dimD;

        //- Force dimensions
        static const dimensionSet dimF;


    // Constructors

        //- Construct from a dictionary and a phase pair
        myTurbulentDispersionModel
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    //- Destructor
    virtual ~myTurbulentDispersionModel();


    // Selectors

        static autoPtr<myTurbulentDispersionModel> New
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    // Member Functions

        //- Turbulent diffusivity
        //  multiplying the gradient of the phase-fraction
        virtual tmp<volScalarField> D() const = 0;    
        virtual tmp<volScalarField> DLam() const;
        virtual tmp<volScalarField> DTurb() const; 
          
        virtual tmp<volVectorField> F() const;
        virtual tmp<volVectorField> FLam() const;
        virtual tmp<volVectorField> FTurb() const; 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
