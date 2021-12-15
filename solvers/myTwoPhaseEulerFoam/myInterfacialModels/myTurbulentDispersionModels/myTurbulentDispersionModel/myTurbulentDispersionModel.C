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

#include "myTurbulentDispersionModel.H"
#include "myPhasePair.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(myTurbulentDispersionModel, 0);
    defineRunTimeSelectionTable(myTurbulentDispersionModel, dictionary);
}

const Foam::dimensionSet Foam::myTurbulentDispersionModel::dimD(1, -1, -2, 0, 0);
const Foam::dimensionSet Foam::myTurbulentDispersionModel::dimF(1, -2, -2, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myTurbulentDispersionModel::myTurbulentDispersionModel
(
    const dictionary& dict,
    const myPhasePair& pair
)
:
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myTurbulentDispersionModel::~myTurbulentDispersionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::myTurbulentDispersionModel::F() const
{
    return D()*fvc::grad(pair_.dispersed());
}



Foam::tmp<Foam::volVectorField>
Foam::myTurbulentDispersionModel::FLam() const
{
    return DLam()*fvc::grad(pair_.dispersed());
}



Foam::tmp<Foam::volVectorField>
Foam::myTurbulentDispersionModel::FTurb() const
{
    return DTurb()*fvc::grad(pair_.dispersed());
}

Foam::tmp<Foam::volScalarField>
Foam::myTurbulentDispersionModel::DLam() const
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
Foam::myTurbulentDispersionModel::DTurb() const
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

// ************************************************************************* //
