/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2017 OpenFOAM Foundation
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

#include "kOmegaSSTSato.H"
#include "fvOptions.H"
#include "myTwoPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
kOmegaSSTSato<BasicMyTurbulenceModel>::kOmegaSSTSato
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
    kOmegaSST<BasicMyTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    gasTurbulencePtr_(nullptr),

    Cmub_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmub",
            this->coeffDict_,
            0.6
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
bool kOmegaSSTSato<BasicMyTurbulenceModel>::read()
{
    if (kOmegaSST<BasicMyTurbulenceModel>::read())
    {
        Cmub_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicMyTurbulenceModel>
const PhaseCompressibleMyTurbulenceModel
<
    typename BasicMyTurbulenceModel::transportModel
>&
kOmegaSSTSato<BasicMyTurbulenceModel>::gasTurbulence() const
{
    if (!gasTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const transportModel& liquid = this->transport();
        const myTwoPhaseSystem& fluid =
            refCast<const myTwoPhaseSystem>(liquid.fluid());
        const transportModel& gas = fluid.otherPhase(liquid);

        gasTurbulencePtr_ =
           &U.db()
           .lookupObject<PhaseCompressibleMyTurbulenceModel<transportModel>>
            (
                IOobject::groupName
                (
                    myTurbulenceModel::propertiesName,
                    gas.name()
                )
            );
    }

    return *gasTurbulencePtr_;
}


template<class BasicMyTurbulenceModel>
void kOmegaSSTSato<BasicMyTurbulenceModel>::correctNut
(
    const volScalarField& S2
)
{
    const PhaseCompressibleMyTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    volScalarField yPlus
    (
        pow(this->betaStar_, 0.25)*this->y_*sqrt(this->k_)/this->nu()
    );

    this->nut_ =
        this->a1_*this->k_
       /max
        (
            this->a1_*this->omega_,
            this->b1_*this->F23()*sqrt(S2)
        )
      + sqr(1 - exp(-yPlus/16.0))
       *Cmub_*gasTurbulence.transport().d()*gasTurbulence.alpha()
       *(mag(this->U_ - gasTurbulence.U()));

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicMyTurbulenceModel::correctNut();
}


template<class BasicMyTurbulenceModel>
void kOmegaSSTSato<BasicMyTurbulenceModel>::correct()
{
    kOmegaSST<BasicMyTurbulenceModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
