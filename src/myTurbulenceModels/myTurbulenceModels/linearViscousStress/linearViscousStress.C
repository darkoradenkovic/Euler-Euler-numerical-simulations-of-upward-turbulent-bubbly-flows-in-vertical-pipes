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

#include "linearViscousStress.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
Foam::linearViscousStress<BasicMyTurbulenceModel>::linearViscousStress
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
    BasicMyTurbulenceModel
    (
        modelName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
bool Foam::linearViscousStress<BasicMyTurbulenceModel>::read()
{
    return BasicMyTurbulenceModel::read();
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::linearViscousStress<BasicMyTurbulenceModel>::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (-(this->alpha_*this->rho_*this->nuEff()))
           *dev(twoSymm(fvc::grad(this->U_)))
        )
    );
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::linearViscousStress<BasicMyTurbulenceModel>::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvc::div((this->alpha_*this->rho_*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*this->rho_*this->nuEff(), U)
    );
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::fvVectorMatrix>
Foam::linearViscousStress<BasicMyTurbulenceModel>::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    return
    (
      - fvc::div((this->alpha_*rho*this->nuEff())*dev2(T(fvc::grad(U))))
      - fvm::laplacian(this->alpha_*rho*this->nuEff(), U)
    );
}


template<class BasicMyTurbulenceModel>
void Foam::linearViscousStress<BasicMyTurbulenceModel>::correct()
{
    BasicMyTurbulenceModel::correct();
}


// ************************************************************************* //
