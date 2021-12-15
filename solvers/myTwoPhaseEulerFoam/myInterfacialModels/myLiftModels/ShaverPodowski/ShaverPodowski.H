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
    Foam::myLiftModels::ShaverPodowski


SourceFiles
    ShaverPodowski.C

\*---------------------------------------------------------------------------*/

#ifndef ShaverPodowski_H
#define ShaverPodowski_H

#include "myLiftModel.H"
#include "myWallDependentModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myPhasePair;

namespace myLiftModels
{

/*---------------------------------------------------------------------------*\
                           Class ShaverPodowski Declaration
\*---------------------------------------------------------------------------*/

class ShaverPodowski
:
    public myLiftModel,
    public myWallDependentModel
{
        const dimensionedScalar CL0_;
    
public:
    
    //- Runtime type information
    TypeName("ShaverPodowski");


    // Constructors

        //- Construct from a dictionary and a phase pair
        ShaverPodowski
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    //- Destructor
    virtual ~ShaverPodowski();


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace myLiftModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
