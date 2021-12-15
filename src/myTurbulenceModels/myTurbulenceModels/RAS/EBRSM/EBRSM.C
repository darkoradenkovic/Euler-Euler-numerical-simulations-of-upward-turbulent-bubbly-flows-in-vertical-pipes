/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010, 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "EBRSM.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
    
template<class BasicMyTurbulenceModel>
void EBRSM<BasicMyTurbulenceModel>::correctNut() 
{

}


template<class BasicMyTurbulenceModel>
tmp<volScalarField> EBRSM<BasicMyTurbulenceModel>::Ls() const
{
    return
        CL_
        *max
        (
            pow(k_,1.5)/epsilon_,
            CEta_*pow(this->nu(),0.75)/pow(epsilon_,0.25)
            
        );
}

template<class BasicMyTurbulenceModel>
tmp<volScalarField> EBRSM<BasicMyTurbulenceModel>::Ts() const
{
    return
        max
        (
            k_/epsilon_,
            CT_*sqrt(this->nu()/epsilon_)
        );
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
EBRSM<BasicMyTurbulenceModel>::EBRSM
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
    ReynoldsStress<RASModel<BasicMyTurbulenceModel>>
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
    
    g1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "g1",
            this->coeffDict_,
            3.4
        )
    ),

    g1Star_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "g1Star",
            this->coeffDict_,
            1.8
        )
    ),

    g3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "g3",
            this->coeffDict_,
            0.8
        )
    ),
    g3Star_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "g3Star",
            this->coeffDict_,
            1.3
        )
    ),
        g4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "g4",
            this->coeffDict_,
            1.25
        )
    ),
        g5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "g5",
            this->coeffDict_,
            0.4
        )
    ),
         
     CMu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CMu",
            this->coeffDict_,
            0.21
        )
    ),
    
    
        sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            1.0
        )
    ),
    
    CT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CT",
            this->coeffDict_,
            6.0
        )
    ),
    
    
    CL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            this->coeffDict_,
            0.133
        )
    ),
    CEta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CEta",
            this->coeffDict_,
            80.0
        )
    ),

    CEps1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CEps1",
            this->coeffDict_,
            1.44
        )
    ),
    
    CEps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CEps2",
            this->coeffDict_,
            1.83
        )
    ),
    
    A1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A1",
            this->coeffDict_,
            0.065
        )
    ),
    
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.15
        )
    ),
    
    
    alphaMH_
    (
        IOobject
        (
            IOobject::groupName("alphaMH", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
        k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*tr(this->R_)
    ),
         epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    alphaMHMin_(dimensionedScalar("alphaMHMin", alphaMHMin_.dimensions(), Zero))
    
{
    if (type == typeName)
    {
        this->printCoeffs(type);

        this->boundNormalStress(this->R_);
        bound(epsilon_, this->epsilonMin_);
        k_ = 0.5*tr(this->R_);
        bound(alphaMH_,alphaMHMin_);
}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
bool EBRSM<BasicMyTurbulenceModel>::read()
{
    if (ReynoldsStress<RASModel<BasicMyTurbulenceModel>>::read())
    {
             g1_.readIfPresent(this->coeffDict());
             g1Star_.readIfPresent(this->coeffDict());
             g3_.readIfPresent(this->coeffDict());
             g3Star_.readIfPresent(this->coeffDict());
             g4_.readIfPresent(this->coeffDict());
             g5_.readIfPresent(this->coeffDict());
             CMu_.readIfPresent(this->coeffDict());
             sigmaK_.readIfPresent(this->coeffDict());
             CT_.readIfPresent(this->coeffDict());
             CL_.readIfPresent(this->coeffDict());
             CEta_.readIfPresent(this->coeffDict());
             CEps1_.readIfPresent(this->coeffDict());
             CEps2_.readIfPresent(this->coeffDict());
             A1_.readIfPresent(this->coeffDict());
             sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


 template<class BasicMyTurbulenceModel>
 tmp<volSymmTensorField> EBRSM<BasicMyTurbulenceModel>::DREff() const
 {
     return tmp<volSymmTensorField>
     (
         new volSymmTensorField
         (
             "DREff",
              (CMu_/sigmaK_*Ts())*this->R_ + I*this->nu()
         )
     );
 }
 
 template<class BasicMyTurbulenceModel>
 tmp<volSymmTensorField> EBRSM<BasicMyTurbulenceModel>::DepsilonEff() const
 {
     return tmp<volSymmTensorField>
     (
         new volSymmTensorField
         (
             "DepsilonEff",
             (CMu_/sigmaEps_*Ts())*this->R_ + I*this->nu()
         )
     );
 }

template<class BasicMyTurbulenceModel>
void EBRSM<BasicMyTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volSymmTensorField& R = this->R_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

  //  ReynoldsStress<RASModel<BasicMyTurbulenceModel>>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField P(-twoSymm(R & gradU));
    volScalarField G(this->GName(), 0.5*mag(tr(P)));
  
    const volScalarField T(this->Ts());
    const volScalarField L(this->Ls());
      
      
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        CEps1_*(1.0 + A1_*(1.0-pow(alphaMH_, 3))*(G/epsilon_))*alpha*rho*G/T
      - fvm::Sp(alpha*rho*CEps2_/T, epsilon_)
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);
             
 
    volSymmTensorField b(dev(R)/(2*k_));
    volSymmTensorField S(symm(gradU));
    volTensorField W(skew(gradU));
    
       
    volVectorField n = fvc::grad(alphaMH_);
    
    n /= mag(n) + dimensionedScalar("nsmall", n.dimensions(), VSMALL);
    
    volSymmTensorField nn = symm(n*n);
    
    tmp<fvSymmTensorMatrix> REqn
    (
       fvm::ddt(alpha, rho, R)
      +fvm::div(alphaRhoPhi, R)
      -fvm::laplacian(alpha*rho*DREff(),R)
      +fvm::Sp(alpha*rho*pow(alphaMH_, 3)*(g1_*epsilon_+g1Star_*G)/(2.0*k_), R)
      +fvm::Sp(alpha*rho*(1.0-pow(alphaMH_, 3))*epsilon_/k_, R)
      ==
       alpha*rho*P
      -((1.0/3.0)*I) * ((pow(alphaMH_, 3)*((2.0-g1_)*epsilon_ - g1Star_*G))*alpha*rho)
      +alpha*rho*((1.0-pow(alphaMH_, 3))*(-5.0*epsilon_/k_*(twoSymm(R & nn) - 0.5 * (R && nn) * (nn + I))))
      +alpha*rho*(pow(alphaMH_, 3)*k_*
      (
          (g3_-g3Star_*mag(b))*dev(S)
          +g4_*(dev(twoSymm(b&S)))
          + g5_*twoSymm(b&W)
      )
      )
      +fvOptions(alpha, rho, R)
    );
   
    REqn.ref().relax();
    fvOptions.constrain(REqn.ref());
    solve(REqn);
    fvOptions.correct(R);

    this->boundNormalStress(R);

    k_ = 0.5*tr(R);
    bound(k_, this->kMin_);

    volScalarField OneOverL2 = 1.0/(L*L);
     
    tmp<fvScalarMatrix> alphaMHEqn
    (
      fvm::Sp(OneOverL2, alphaMH_)-fvm::laplacian(alphaMH_) == OneOverL2
    );
    
    
    alphaMHEqn.ref().relax();
    fvOptions.constrain(alphaMHEqn.ref());
    solve(alphaMHEqn);   
    fvOptions.correct(alphaMH_); 
    bound(alphaMH_, alphaMHMin_); 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
