/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2015-2017 OpenFOAM Foundation
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

#include "ThermalDiffusivity.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
Foam::ThermalDiffusivity<BasicMyTurbulenceModel>::ThermalDiffusivity
(
    const word& type,
    const alphaField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicMyTurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
Foam::autoPtr<Foam::ThermalDiffusivity<BasicMyTurbulenceModel>>
Foam::ThermalDiffusivity<BasicMyTurbulenceModel>::New
(
    const alphaField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    return autoPtr<ThermalDiffusivity>
    (
        static_cast<ThermalDiffusivity*>(
        BasicMyTurbulenceModel::New
        (
            alpha,
            rho,
            U,
            alphaRhoPhi,
            phi,
            transport,
            propertiesName
        ).ptr())
    );
}


template<class BasicMyTurbulenceModel>
Foam::autoPtr<Foam::ThermalDiffusivity<BasicMyTurbulenceModel>>
Foam::ThermalDiffusivity<BasicMyTurbulenceModel>::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    return autoPtr<ThermalDiffusivity>
    (
        static_cast<ThermalDiffusivity*>(
        BasicMyTurbulenceModel::New
        (
            rho,
            U,
            phi,
            transport,
            propertiesName
        ).ptr())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::ThermalDiffusivity<BasicMyTurbulenceModel>::alphat() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("alphat", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimDensity*dimViscosity, Zero)
    );
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::ThermalDiffusivity<BasicMyTurbulenceModel>::alphat
(
    const label patchi
) const
{
    return tmp<scalarField>
    (
        new scalarField(this->mesh_.boundary()[patchi].size(), Zero)
    );
}


// ************************************************************************* //
