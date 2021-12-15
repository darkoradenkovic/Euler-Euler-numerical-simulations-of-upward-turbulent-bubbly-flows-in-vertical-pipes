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

\*---------------------------------------------------------------------------*/

#include "TomiyamaAspectRatio.H"
#include "myPhasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myAspectRatioModels
{
    defineTypeNameAndDebug(TomiyamaAspectRatio, 0);
    addToRunTimeSelectionTable
    (
        myAspectRatioModel,
        TomiyamaAspectRatio,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myAspectRatioModels::TomiyamaAspectRatio::TomiyamaAspectRatio
(
    const dictionary& dict,
    const myPhasePair& pair
)
:
    VakhrushevEfremov(dict, pair),
    myWallDependentModel(pair.phase1().mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myAspectRatioModels::TomiyamaAspectRatio::~TomiyamaAspectRatio()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::myAspectRatioModels::TomiyamaAspectRatio::E() const
{
    return
        VakhrushevEfremov::E()
       *max
       (
           scalar(1) - 0.35*yWall()/pair_.dispersed().d(),
           scalar(0.65)
       );
}


// ************************************************************************* //
