/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

#include "continuousGasKEpsilon.H"
#include "fvOptions.H"
#include "myTwoPhaseSystem.H"
#include "myVirtualMassModel.H"
#include "myDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
continuousGasKEpsilon<BasicMyTurbulenceModel>::continuousGasKEpsilon
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
    kEpsilon<BasicMyTurbulenceModel>
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

    liquidTurbulencePtr_(nullptr),

    nutEff_
    (
        IOobject
        (
            IOobject::groupName("nutEff", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->nut_
    ),

    alphaInversion_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaInversion",
            this->coeffDict_,
            0.7
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
bool continuousGasKEpsilon<BasicMyTurbulenceModel>::read()
{
    if (kEpsilon<BasicMyTurbulenceModel>::read())
    {
        alphaInversion_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicMyTurbulenceModel>
void continuousGasKEpsilon<BasicMyTurbulenceModel>::correctNut()
{
    kEpsilon<BasicMyTurbulenceModel>::correctNut();

    const myTurbulenceModel& liquidTurbulence = this->liquidTurbulence();
    const transportModel& gas = this->transport();
    const myTwoPhaseSystem& fluid = refCast<const myTwoPhaseSystem>(gas.fluid());
    const transportModel& liquid = fluid.otherPhase(gas);

     const myVirtualMassModel& virtualMass =
        fluid.lookupSubModel<myVirtualMassModel>(gas, liquid);

    volScalarField thetal(liquidTurbulence.k()/liquidTurbulence.epsilon());
    volScalarField rhodv(gas.rho() + virtualMass.Cvm()*liquid.rho());
    volScalarField thetag((rhodv/(18*liquid.rho()*liquid.nu()))*sqr(gas.d()));
    volScalarField expThetar
    (
        min
        (
            exp(min(thetal/thetag, scalar(50))),
            scalar(1)
        )
    );
    volScalarField omega((1 - expThetar)/(1 + expThetar));

    nutEff_ = omega*liquidTurbulence.nut();
    fv::options::New(this->mesh_).correct(nutEff_);
}


template<class BasicMyTurbulenceModel>
const myTurbulenceModel&
continuousGasKEpsilon<BasicMyTurbulenceModel>::liquidTurbulence() const
{
    if (!liquidTurbulencePtr_)
    {
        const volVectorField& U = this->U_;

        const transportModel& gas = this->transport();
        const myTwoPhaseSystem& fluid =
            refCast<const myTwoPhaseSystem>(gas.fluid());
        const transportModel& liquid = fluid.otherPhase(gas);

        liquidTurbulencePtr_ =
           &U.db().lookupObject<myTurbulenceModel>
            (
                IOobject::groupName
                (
                    myTurbulenceModel::propertiesName,
                    liquid.name()
                )
            );
    }

    return *liquidTurbulencePtr_;
}


template<class BasicMyTurbulenceModel>
tmp<Foam::volScalarField>
continuousGasKEpsilon<BasicMyTurbulenceModel>::nuEff() const
{
    volScalarField blend
    (
        max
        (
            min
            (
                (this->alpha_ - scalar(0.5))/(alphaInversion_ - 0.5),
                scalar(1)
            ),
            scalar(0)
        )
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("nuEff", this->alphaRhoPhi_.group()),
            blend*this->nut_
          + (1.0 - blend)*rhoEff()*nutEff_/this->transport().rho()
          + this->nu()
        )
    );
}


template<class BasicMyTurbulenceModel>
tmp<Foam::volScalarField>
continuousGasKEpsilon<BasicMyTurbulenceModel>::rhoEff() const
{
    const transportModel& gas = this->transport();
    const myTwoPhaseSystem& fluid = refCast<const myTwoPhaseSystem>(gas.fluid());
    const transportModel& liquid = fluid.otherPhase(gas);

    const myVirtualMassModel& virtualMass =
        fluid.lookupSubModel<myVirtualMassModel>(gas, liquid);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("rhoEff", this->alphaRhoPhi_.group()),
            gas.rho() + (virtualMass.Cvm() + 3.0/20.0)*liquid.rho()
        )
    );
}


template<class BasicMyTurbulenceModel>
tmp<volScalarField>
continuousGasKEpsilon<BasicMyTurbulenceModel>::phaseTransferCoeff() const
{
    const volVectorField& U = this->U_;
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;

    const myTurbulenceModel& liquidTurbulence = this->liquidTurbulence();

    return
    (
        max(alphaInversion_ - alpha, scalar(0))
       *rho
       *min
        (
            liquidTurbulence.epsilon()/liquidTurbulence.k(),
            1.0/U.time().deltaT()
        )
    );
}


template<class BasicMyTurbulenceModel>
tmp<fvScalarMatrix>
continuousGasKEpsilon<BasicMyTurbulenceModel>::kSource() const
{
    const myTurbulenceModel& liquidTurbulence = this->liquidTurbulence();
    const volScalarField phaseTransferCoeff(this->phaseTransferCoeff());

    return
        phaseTransferCoeff*liquidTurbulence.k()
      - fvm::Sp(phaseTransferCoeff, this->k_);
}


template<class BasicMyTurbulenceModel>
tmp<fvScalarMatrix>
continuousGasKEpsilon<BasicMyTurbulenceModel>::epsilonSource() const
{
    const myTurbulenceModel& liquidTurbulence = this->liquidTurbulence();
    const volScalarField phaseTransferCoeff(this->phaseTransferCoeff());

    return
        phaseTransferCoeff*liquidTurbulence.epsilon()
      - fvm::Sp(phaseTransferCoeff, this->epsilon_);
}


template<class BasicMyTurbulenceModel>
tmp<volSymmTensorField>
continuousGasKEpsilon<BasicMyTurbulenceModel>::R() const
{
    tmp<volScalarField> tk(this->k());

    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", this->alphaRhoPhi_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*tk() - (nutEff_)*dev(twoSymm(fvc::grad(this->U_))),
            tk().boundaryField().types()
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
