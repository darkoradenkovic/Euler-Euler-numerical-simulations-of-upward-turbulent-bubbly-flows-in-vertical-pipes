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

#include "nonlinearEddyViscosity.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
Foam::nonlinearEddyViscosity<BasicMyTurbulenceModel>::nonlinearEddyViscosity
(
    const word& modelName,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    eddyViscosity<BasicMyTurbulenceModel>
    (
        modelName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    nonlinearStress_
    (
        IOobject
        (
            IOobject::groupName("nonlinearStress", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedSymmTensor(sqr(dimVelocity), Zero)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::nonlinearEddyViscosity<BasicMyTurbulenceModel>::R() const
{
    tmp<volSymmTensorField> tR
    (
        eddyViscosity<BasicMyTurbulenceModel>::R()
    );
    tR.ref() += nonlinearStress_;
    return tR;
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::nonlinearEddyViscosity<BasicMyTurbulenceModel>::devRhoReff() const
{
    tmp<volSymmTensorField> tdevRhoReff
    (
        eddyViscosity<BasicMyTurbulenceModel>::devRhoReff()
    );
    tdevRhoReff.ref() += this->rho_*nonlinearStress_;
    return tdevRhoReff;
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::nonlinearEddyViscosity<BasicMyTurbulenceModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
        fvc::div(this->rho_*nonlinearStress_)
      + eddyViscosity<BasicMyTurbulenceModel>::divDevRhoReff(U)
    );
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::nonlinearEddyViscosity<BasicMyTurbulenceModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return
    (
        fvc::div(rho*nonlinearStress_)
      + eddyViscosity<BasicMyTurbulenceModel>::divDevRhoReff(rho, U)
    );
}


// ************************************************************************* //
