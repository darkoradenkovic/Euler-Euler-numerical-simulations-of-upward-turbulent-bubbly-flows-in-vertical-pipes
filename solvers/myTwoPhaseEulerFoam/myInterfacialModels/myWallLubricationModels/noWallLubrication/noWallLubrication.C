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

#include "noWallLubrication.H"
#include "myPhasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myWallLubricationModels
{
    defineTypeNameAndDebug(noWallLubrication, 0);
    addToRunTimeSelectionTable
    (
        myWallLubricationModel,
        noWallLubrication,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myWallLubricationModels::noWallLubrication::noWallLubrication
(
    const dictionary& dict,
    const myPhasePair& pair
)
:
    myWallLubricationModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myWallLubricationModels::noWallLubrication::~noWallLubrication()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::myWallLubricationModels::noWallLubrication::Fi() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return tmp<volVectorField>::New
    (
        IOobject
        (
            "noWallLubrication:Fi",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedVector(dimF, Zero)
    );
}


Foam::tmp<Foam::volVectorField>
Foam::myWallLubricationModels::noWallLubrication::F() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return tmp<volVectorField>::New
    (
        IOobject
        (
            "noWallLubrication:F",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedVector(dimF, Zero)
    );
}


// ************************************************************************* //
