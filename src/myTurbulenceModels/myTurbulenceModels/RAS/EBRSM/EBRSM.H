/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

Class
    Foam::RASModels::EBRSM

Group
    grpRASTurbulence

Description
    Launder, Reece and Rodi Reynolds-stress turbulence model for
    incompressible and compressible flows.

    Reference:
    \verbatim
        Launder, B. E., Reece, G. J., & Rodi, W. (1975).
        Progress in the development of a Reynolds-stress turbulence closure.
        Journal of fluid mechanics, 68(03), 537-566.
    \endverbatim

    Including the recommended generalized gradient diffusion model of
    Daly and Harlow:
    \verbatim
        Daly, B. J., & Harlow, F. H. (1970).
        Transport equations in turbulence.
        Physics of Fluids (1958-1988), 13(11), 2634-2649.
    \endverbatim

    Optional Gibson-Launder wall-reflection is also provided:
    \verbatim
        Gibson, M. M., & Launder, B. E. (1978).
        Ground effects on pressure fluctuations in the
        atmospheric boundary layer.
        Journal of Fluid Mechanics, 86(03), 491-511.
    \endverbatim

    The default model coefficients are:
    \verbatim
        EBRSMCoeffs
        {
            Cmu             0.09;
            C1              1.8;
            C2              0.6;
            Ceps1           1.44;
            Ceps2           1.92;
            Cs              0.25;
            Ceps            0.15;

            wallReflection  yes;
            kappa           0.41
            Cref1           0.5;
            Cref2           0.3;

            couplingFactor  0.0;
        }
    \endverbatim

SourceFiles
    EBRSM.C

\*---------------------------------------------------------------------------*/

#ifndef EBRSM_H
#define EBRSM_H

#include "RASModel.H"
#include "ReynoldsStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class EBRSM Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMyTurbulenceModel>
class EBRSM
:
    public ReynoldsStress<RASModel<BasicMyTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        EBRSM(const EBRSM&) = delete;

        //- No copy assignment
        void operator=(const EBRSM&) = delete;


protected:

    // Protected data

        // Model coefficients
            dimensionedScalar g1_;
            dimensionedScalar g1Star_;
       //     dimensionedScalar g2_;
            dimensionedScalar g3_;
            dimensionedScalar g3Star_;
            dimensionedScalar g4_;
            dimensionedScalar g5_;
            dimensionedScalar CMu_;
            dimensionedScalar sigmaK_;
            dimensionedScalar CT_;
            dimensionedScalar CL_;
            dimensionedScalar CEta_;
            dimensionedScalar CEps1_;
            dimensionedScalar CEps2_;
            dimensionedScalar A1_;
            dimensionedScalar sigmaEps_;     
            
         

        // Fields
            volScalarField alphaMH_;
            volScalarField k_;
            volScalarField epsilon_;
            dimensionedScalar alphaMHMin_;
    // Protected Member Functions

        //- Update the eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicMyTurbulenceModel::alphaField alphaField;
    typedef typename BasicMyTurbulenceModel::rhoField rhoField;
    typedef typename BasicMyTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("EBRSM");


    // Constructors

        //- Construct from components
        EBRSM
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = myTurbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~EBRSM()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }
        
//                 //- Return the effective stress tensor
//         virtual tmp<volSymmTensorField> devRhoReff() const
//         {
// //           fvc::laplacian
// //           (
// //                  this->alpha_*rho*this->nut(),
// //                  U,
// //                 "laplacian(nuEff,U)"
// //           )
//             return
//             (
//             fvc::div(this->alpha_*this->rho_*this->R_)
// //          - fvc::div(this->alpha_*rho*this->nu()*dev2(T(fvc::grad(U))))
//           - fvm::laplacian(this->alpha_*this->rho_*this->nu(), this->U_)
//             );
//         }
        
        

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the effective diffusivity for R
        tmp<volSymmTensorField> DREff() const;

        //- Return the effective diffusivity for epsilon
        tmp<volSymmTensorField> DepsilonEff() const;

        //- Solve the turbulence equations and correct eddy-Viscosity and
        //  related properties
                //- Return time scale, Ts
        tmp<volScalarField> Ts() const;

        //- Return length scale, Ls
        tmp<volScalarField> Ls() const;
        
        
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EBRSM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
