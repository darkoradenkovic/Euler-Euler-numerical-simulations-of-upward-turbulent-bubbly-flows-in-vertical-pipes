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

#include "ShaverPodowski.H"
#include "myPhasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myLiftModels
{
    defineTypeNameAndDebug(ShaverPodowski, 0);
    addToRunTimeSelectionTable(myLiftModel, ShaverPodowski, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myLiftModels::ShaverPodowski::ShaverPodowski
(
    const dictionary& dict,
    const myPhasePair& pair
)
:
    myLiftModel(dict, pair),
    myWallDependentModel(pair.phase1().mesh()),
    CL0_("CL0", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myLiftModels::ShaverPodowski::~ShaverPodowski()
{}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::myLiftModels::ShaverPodowski::Cl() const
{
   
    const volScalarField& y(yWall());
    volScalarField yTilde(y/(pair_.dispersed().d()));

    return
        neg(yTilde - 0.5)*0
      + pos0(yTilde - 0.5)*pos0(1.0-yTilde)*CL0_*(3.0*sqr(2.0*yTilde-1.0) - 2.0*pow((2.0*yTilde-1.0),3))
      + neg(1.0-yTilde)*CL0_;
}

// ************************************************************************* //
