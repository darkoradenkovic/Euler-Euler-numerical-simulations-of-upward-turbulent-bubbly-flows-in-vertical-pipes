/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013-2014 OpenFOAM Foundation
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

#include "myIncompressibleMyTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(myIncompressibleMyTurbulenceModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myIncompressibleMyTurbulenceModel::myIncompressibleMyTurbulenceModel
(
    const geometricOneField&,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const word& propertiesName
)
:
    myTurbulenceModel
    (
        U,
        alphaRhoPhi,
        phi,
        propertiesName
    )
{}


Foam::tmp<Foam::volScalarField>
Foam::myIncompressibleMyTurbulenceModel::mu() const
{
    return nu();
}


Foam::tmp<Foam::scalarField>
Foam::myIncompressibleMyTurbulenceModel::mu(const label patchi) const
{
    return nu(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::myIncompressibleMyTurbulenceModel::mut() const
{
    return nut();
}


Foam::tmp<Foam::scalarField>
Foam::myIncompressibleMyTurbulenceModel::mut(const label patchi) const
{
    return nut(patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::myIncompressibleMyTurbulenceModel::muEff() const
{
    return nuEff();
}


Foam::tmp<Foam::scalarField>
Foam::myIncompressibleMyTurbulenceModel::muEff(const label patchi) const
{
    return nuEff(patchi);
}


// ************************************************************************* //
