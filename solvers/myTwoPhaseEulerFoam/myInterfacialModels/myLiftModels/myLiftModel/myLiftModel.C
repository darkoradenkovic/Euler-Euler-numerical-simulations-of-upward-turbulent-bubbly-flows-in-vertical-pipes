/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "myLiftModel.H"
#include "myPhasePair.H"
#include "fvcCurl.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(myLiftModel, 0);
    defineRunTimeSelectionTable(myLiftModel, dictionary);
}

const Foam::dimensionSet Foam::myLiftModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myLiftModel::myLiftModel
(
    const dictionary& dict,
    const myPhasePair& pair
)
:
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myLiftModel::~myLiftModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField> Foam::myLiftModel::Fi() const
{
    return
        Cl()
       *pair_.continuous().rho()
       *(
            pair_.Ur() ^ fvc::curl(pair_.continuous().U())
        );
}


Foam::tmp<Foam::volVectorField> Foam::myLiftModel::F() const
{
    return pair_.dispersed()*Fi();
}


Foam::tmp<Foam::surfaceScalarField> Foam::myLiftModel::Ff() const
{
    return fvc::interpolate(pair_.dispersed())*fvc::flux(Fi());
}


// ************************************************************************* //
