/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

#include "NicenoKEqn.H"
#include "fvOptions.H"
#include "myTwoPhaseSystem.H"
#include "myDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
NicenoKEqn<BasicMyTurbulenceModel>::NicenoKEqn
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
    kEqn<BasicMyTurbulenceModel>
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

    alphaInversion_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaInversion",
            this->coeffDict_,
            0.3
        )
    ),

    Cp_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cp",
            this->coeffDict_,
            this->Ck_.value()
        )
    ),

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
bool NicenoKEqn<BasicMyTurbulenceModel>::read()
{
    if (kEqn<BasicMyTurbulenceModel>::read())
    {
        alphaInversion_.readIfPresent(this->coeffDict());
        Cp_.readIfPresent(this->coeffDict());
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
NicenoKEqn<BasicMyTurbulenceModel>::gasTurbulence() const
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
void NicenoKEqn<BasicMyTurbulenceModel>::correctNut()
{
    const PhaseCompressibleMyTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    this->nut_ =
        this->Ck_*sqrt(this->k_)*this->delta()
      + Cmub_*gasTurbulence.transport().d()*gasTurbulence.alpha()
       *(mag(this->U_ - gasTurbulence.U()));

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicMyTurbulenceModel::correctNut();
}


template<class BasicMyTurbulenceModel>
tmp<volScalarField> NicenoKEqn<BasicMyTurbulenceModel>::bubbleG() const
{
    const PhaseCompressibleMyTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    const transportModel& liquid = this->transport();
    const myTwoPhaseSystem& fluid =
        refCast<const myTwoPhaseSystem>(liquid.fluid());
    const transportModel& gas = fluid.otherPhase(liquid);

    const myDragModel& drag = fluid.lookupSubModel<myDragModel>(gas, liquid);

    volScalarField magUr(mag(this->U_ - gasTurbulence.U()));

    tmp<volScalarField> bubbleG
    (
        Cp_*sqr(magUr)*drag.K()/liquid.rho()
    );

    return bubbleG;
}


template<class BasicMyTurbulenceModel>
tmp<volScalarField>
NicenoKEqn<BasicMyTurbulenceModel>::phaseTransferCoeff() const
{
    const volVectorField& U = this->U_;
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const myTurbulenceModel& gasTurbulence = this->gasTurbulence();

    return
    (
        max(alphaInversion_ - alpha, scalar(0))
       *rho
       *min
        (
            this->Ce_*sqrt(gasTurbulence.k())/this->delta(),
            1.0/U.time().deltaT()
        )
    );
}


template<class BasicMyTurbulenceModel>
tmp<fvScalarMatrix> NicenoKEqn<BasicMyTurbulenceModel>::kSource() const
{
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const PhaseCompressibleMyTurbulenceModel<transportModel>& gasTurbulence =
        this->gasTurbulence();

    const volScalarField phaseTransferCoeff(this->phaseTransferCoeff());

    return
        alpha*rho*bubbleG()
      + phaseTransferCoeff*gasTurbulence.k()
      - fvm::Sp(phaseTransferCoeff, this->k_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
