/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2015 OpenFOAM Foundation
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

#include "kOmegaSSTDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
void kOmegaSSTDES<BasicMyTurbulenceModel>::correctNut(const volScalarField& S2)
{
    // Correct the turbulence viscosity
    kOmegaSSTBase<DESModel<BasicMyTurbulenceModel>>::correctNut(S2);

    // Correct the turbulence thermal diffusivity
    BasicMyTurbulenceModel::correctNut();
}


template<class BasicMyTurbulenceModel>
void kOmegaSSTDES<BasicMyTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicMyTurbulenceModel>
tmp<volScalarField> kOmegaSSTDES<BasicMyTurbulenceModel>::dTilda
(
    const volScalarField& magGradU,
    const volScalarField& CDES
) const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;

    return min(CDES*this->delta(), sqrt(k)/(this->betaStar_*omega));
}


template<class BasicMyTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTDES<BasicMyTurbulenceModel>::epsilonByk
(
    const volScalarField& F1,
    const volTensorField& gradU
) const
{
    volScalarField CDES(this->CDES(F1));
    return sqrt(this->k_())/dTilda(mag(gradU), CDES)()();
}


template<class BasicMyTurbulenceModel>
tmp<volScalarField::Internal> kOmegaSSTDES<BasicMyTurbulenceModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    return GbyNu0; // Unlimited
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
kOmegaSSTDES<BasicMyTurbulenceModel>::kOmegaSSTDES
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSSTBase<DESModel<BasicMyTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    CDESkom_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDESkom",
            this->coeffDict_,
            0.82
        )
    ),
    CDESkeps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDESkeps",
            this->coeffDict_,
            0.60
        )
    )
{
    correctNut();

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
bool kOmegaSSTDES<BasicMyTurbulenceModel>::read()
{
    if (kOmegaSSTBase<DESModel<BasicMyTurbulenceModel>>::read())
    {
        kappa_.readIfPresent(this->coeffDict());
        CDESkom_.readIfPresent(this->coeffDict());
        CDESkeps_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicMyTurbulenceModel>
tmp<volScalarField> kOmegaSSTDES<BasicMyTurbulenceModel>::LESRegion() const
{
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const volVectorField& U = this->U_;

    const volScalarField CDkOmega
    (
        (2*this->alphaOmega2_)*(fvc::grad(k) & fvc::grad(omega))/omega
    );

    const volScalarField F1(this->F1(CDkOmega));

    tmp<volScalarField> tLESRegion
    (
        new volScalarField
        (
            IOobject
            (
                "DES::LESRegion",
                this->mesh_.time().timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            neg
            (
                dTilda
                (
                    mag(fvc::grad(U)),
                    F1*CDESkom_ + (1 - F1)*CDESkeps_
                )
              - sqrt(k)/(this->betaStar_*omega)
            )
        )
    );

    return tLESRegion;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
