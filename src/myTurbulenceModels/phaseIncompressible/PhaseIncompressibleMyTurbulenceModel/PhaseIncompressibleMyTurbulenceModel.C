/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "PhaseIncompressibleMyTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::PhaseIncompressibleMyTurbulenceModel<TransportModel>::
PhaseIncompressibleMyTurbulenceModel
(
    const word& type,
    const volScalarField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transportModel,
    const word& propertiesName
)
:
    MyTurbulenceModel
    <
        volScalarField,
        geometricOneField,
        myIncompressibleMyTurbulenceModel,
        TransportModel
    >
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transportModel,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::PhaseIncompressibleMyTurbulenceModel<TransportModel>>
Foam::PhaseIncompressibleMyTurbulenceModel<TransportModel>::New
(
    const volScalarField& alpha,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transportModel,
    const word& propertiesName
)
{
    return autoPtr<PhaseIncompressibleMyTurbulenceModel>
    (
        static_cast<PhaseIncompressibleMyTurbulenceModel*>(
        MyTurbulenceModel
        <
            volScalarField,
            geometricOneField,
            myIncompressibleMyTurbulenceModel,
            TransportModel
        >::New
        (
            alpha,
            geometricOneField(),
            U,
            alphaRhoPhi,
            phi,
            transportModel,
            propertiesName
        ).ptr())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volScalarField>
Foam::PhaseIncompressibleMyTurbulenceModel<TransportModel>::pPrime() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("pPrime", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimPressure, Zero)
    );
}


template<class TransportModel>
Foam::tmp<Foam::surfaceScalarField>
Foam::PhaseIncompressibleMyTurbulenceModel<TransportModel>::pPrimef() const
{
    return tmp<surfaceScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("pPrimef", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar(dimPressure, Zero)
    );
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::PhaseIncompressibleMyTurbulenceModel<TransportModel>::devReff() const
{
    return devRhoReff();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PhaseIncompressibleMyTurbulenceModel<TransportModel>::divDevReff
(
    volVectorField& U
) const
{
    return divDevRhoReff(U);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::PhaseIncompressibleMyTurbulenceModel<TransportModel>::devRhoReff() const
{
    NotImplemented;

    return devReff();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::PhaseIncompressibleMyTurbulenceModel<TransportModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    NotImplemented;

    return divDevReff(U);
}


// ************************************************************************* //
