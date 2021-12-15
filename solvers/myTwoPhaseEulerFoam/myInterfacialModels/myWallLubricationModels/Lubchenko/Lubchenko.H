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
    Foam::myWallLubricationModels::Lubchenko

Description
    Wall lubrication model of Lubchenko.

    \endverbatim

SourceFiles
    Lubchenko.C

\*---------------------------------------------------------------------------*/

#ifndef Lubchenko_H
#define Lubchenko_H

#include "myWallLubricationModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myPhasePair;

namespace myWallLubricationModels
{

/*---------------------------------------------------------------------------*\
                            Class Lubchenko Declaration
\*---------------------------------------------------------------------------*/

class Lubchenko
:
    public myWallLubricationModel
{
    // Private data

        //- Schmidt number
        const dimensionedScalar sigma_;

        //- Residual phase fraction
        const dimensionedScalar residualAlpha_;

public:

    //- Runtime type information
    TypeName("Lubchenko");


    // Constructors

        //- Construct from components
        Lubchenko
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    //- Destructor
    virtual ~Lubchenko();


    // Member Functions

        //- Return phase-intensive wall lubrication force
        tmp<volVectorField> Fi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace myWallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
