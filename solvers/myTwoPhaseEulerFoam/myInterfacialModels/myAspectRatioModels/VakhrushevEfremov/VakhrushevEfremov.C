/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "VakhrushevEfremov.H"
#include "myPhasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myAspectRatioModels
{
    defineTypeNameAndDebug(VakhrushevEfremov, 0);
    addToRunTimeSelectionTable
    (
        myAspectRatioModel,
        VakhrushevEfremov,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myAspectRatioModels::VakhrushevEfremov::VakhrushevEfremov
(
    const dictionary& dict,
    const myPhasePair& pair
)
:
    myAspectRatioModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myAspectRatioModels::VakhrushevEfremov::~VakhrushevEfremov()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::myAspectRatioModels::VakhrushevEfremov::E() const
{
    volScalarField Ta(pair_.Ta());

    return
        neg(Ta - scalar(1))*scalar(1)
      + pos0(Ta - scalar(1))*neg(Ta - scalar(39.8))
       *pow3(0.81 + 0.206*tanh(1.6 - 2*log10(max(Ta, scalar(1)))))
      + pos0(Ta - scalar(39.8))*0.24;
}


// ************************************************************************* //
