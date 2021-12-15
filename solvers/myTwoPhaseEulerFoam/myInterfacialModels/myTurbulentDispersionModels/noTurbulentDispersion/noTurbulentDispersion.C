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

#include "noTurbulentDispersion.H"
#include "myPhasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myTurbulentDispersionModels
{
    defineTypeNameAndDebug(noTurbulentDispersion, 0);
    addToRunTimeSelectionTable
    (
        myTurbulentDispersionModel,
        noTurbulentDispersion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myTurbulentDispersionModels::noTurbulentDispersion::noTurbulentDispersion
(
    const dictionary& dict,
    const myPhasePair& pair
)
:
    myTurbulentDispersionModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myTurbulentDispersionModels::noTurbulentDispersion::
~noTurbulentDispersion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::myTurbulentDispersionModels::noTurbulentDispersion::D() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "zero",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimD, Zero)
    );
}

Foam::tmp<Foam::volScalarField>
Foam::myTurbulentDispersionModels::noTurbulentDispersion::DLam() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "zero",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimD, Zero)
    );
}

Foam::tmp<Foam::volVectorField>
Foam::myTurbulentDispersionModels::noTurbulentDispersion::F() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return tmp<volVectorField>::New
    (
        IOobject
        (
            "zero",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector(dimF, Zero)
    );
}


// ************************************************************************* //
