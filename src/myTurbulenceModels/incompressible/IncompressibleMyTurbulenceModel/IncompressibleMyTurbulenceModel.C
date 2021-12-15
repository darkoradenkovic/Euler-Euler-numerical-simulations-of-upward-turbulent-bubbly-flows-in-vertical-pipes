/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "IncompressibleMyTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::IncompressibleMyTurbulenceModel<TransportModel>::
IncompressibleMyTurbulenceModel
(
    const word& type,
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const TransportModel& transport,
    const word& propertiesName
)
:
    MyTurbulenceModel
    <
        geometricOneField,
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
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class TransportModel>
Foam::autoPtr<Foam::IncompressibleMyTurbulenceModel<TransportModel>>
Foam::IncompressibleMyTurbulenceModel<TransportModel>::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const TransportModel& transport,
    const word& propertiesName
)
{
    return autoPtr<IncompressibleMyTurbulenceModel>
    (
        static_cast<IncompressibleMyTurbulenceModel*>(
        MyTurbulenceModel
        <
            geometricOneField,
            geometricOneField,
            myIncompressibleMyTurbulenceModel,
            TransportModel
        >::New
        (
            geometricOneField(),
            geometricOneField(),
            U,
            phi,
            phi,
            transport,
            propertiesName
        ).ptr())
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleMyTurbulenceModel<TransportModel>::devReff() const
{
    return devRhoReff();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleMyTurbulenceModel<TransportModel>::divDevReff
(
    volVectorField& U
) const
{
    return divDevRhoReff(U);
}


template<class TransportModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::IncompressibleMyTurbulenceModel<TransportModel>::devRhoReff() const
{
    NotImplemented;

    return devReff();
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleMyTurbulenceModel<TransportModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    NotImplemented;

    return divDevReff(U);
}


template<class TransportModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::IncompressibleMyTurbulenceModel<TransportModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    NotImplemented;

    return divDevReff(U);
}


// ************************************************************************* //
