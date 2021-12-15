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

#include "Stokes.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
Stokes<BasicMyTurbulenceModel>::Stokes
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    linearViscousStress<laminarModel<BasicMyTurbulenceModel>>
    (
        typeName,
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
const dictionary&
Stokes<BasicMyTurbulenceModel>::coeffDict() const
{
    return dictionary::null;
}


template<class BasicMyTurbulenceModel>
bool Stokes<BasicMyTurbulenceModel>::read()
{
    return true;
}


template<class BasicMyTurbulenceModel>
tmp<volScalarField>
Stokes<BasicMyTurbulenceModel>::nut() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("nut", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(dimViscosity, Zero)
    );
}


template<class BasicMyTurbulenceModel>
tmp<scalarField>
Stokes<BasicMyTurbulenceModel>::nut
(
    const label patchi
) const
{
    return tmp<scalarField>
    (
        new scalarField(this->mesh_.boundary()[patchi].size(), Zero)
    );
}


template<class BasicMyTurbulenceModel>
tmp<volScalarField>
Stokes<BasicMyTurbulenceModel>::nuEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("nuEff", this->alphaRhoPhi_.group()), this->nu()
        )
    );
}


template<class BasicMyTurbulenceModel>
tmp<scalarField>
Stokes<BasicMyTurbulenceModel>::nuEff
(
    const label patchi
) const
{
    return this->nu(patchi);
}


template<class BasicMyTurbulenceModel>
tmp<volScalarField>
Stokes<BasicMyTurbulenceModel>::k() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("k", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(sqr(this->U_.dimensions()), Zero)
    );
}


template<class BasicMyTurbulenceModel>
tmp<volScalarField>
Stokes<BasicMyTurbulenceModel>::epsilon() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(sqr(this->U_.dimensions())/dimTime, Zero)
    );
}


template<class BasicMyTurbulenceModel>
tmp<volSymmTensorField>
Stokes<BasicMyTurbulenceModel>::R() const
{
    return tmp<volSymmTensorField>::New
    (
        IOobject
        (
            IOobject::groupName("R", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedSymmTensor(sqr(this->U_.dimensions()), Zero)
    );
}


template<class BasicMyTurbulenceModel>
void Stokes<BasicMyTurbulenceModel>::correct()
{
    laminarModel<BasicMyTurbulenceModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarModels
} // End namespace Foam

// ************************************************************************* //
