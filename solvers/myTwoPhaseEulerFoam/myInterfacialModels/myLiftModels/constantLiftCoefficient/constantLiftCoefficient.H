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
    Foam::myLiftModels::constantLiftCoefficient

Description
    Constant coefficient lift model.

SourceFiles
    constantLiftCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef constantLiftCoefficient_H
#define constantLiftCoefficient_H

#include "myLiftModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myPhasePair;

namespace myLiftModels
{

/*---------------------------------------------------------------------------*\
                           Class constantLiftCoefficient Declaration
\*---------------------------------------------------------------------------*/

class constantLiftCoefficient
:
    public myLiftModel
{
    // Private data

        //- Constant lift coefficient
        const dimensionedScalar Cl_;


public:

    //- Runtime type information
    TypeName("constantCoefficient");


    // Constructors

        //- Construct from a dictionary and a phase pair
        constantLiftCoefficient
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    //- Destructor
    virtual ~constantLiftCoefficient();


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
