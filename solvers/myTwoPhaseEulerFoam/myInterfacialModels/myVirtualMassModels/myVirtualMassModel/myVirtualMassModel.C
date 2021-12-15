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

#include "myVirtualMassModel.H"
#include "myPhasePair.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(myVirtualMassModel, 0);
    defineRunTimeSelectionTable(myVirtualMassModel, dictionary);
}

const Foam::dimensionSet Foam::myVirtualMassModel::dimK(dimDensity);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myVirtualMassModel::myVirtualMassModel
(
    const dictionary& dict,
    const myPhasePair& pair,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    pair_(pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myVirtualMassModel::~myVirtualMassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::myVirtualMassModel::Ki() const
{
    return Cvm()*pair_.continuous().rho();
}


Foam::tmp<Foam::volScalarField> Foam::myVirtualMassModel::K() const
{
    return pair_.dispersed()*Ki();
}


Foam::tmp<Foam::surfaceScalarField> Foam::myVirtualMassModel::Kf() const
{
    return
        fvc::interpolate(pair_.dispersed())*fvc::interpolate(Ki());
}


bool Foam::myVirtualMassModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
