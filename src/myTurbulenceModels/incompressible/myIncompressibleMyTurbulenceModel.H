/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013-2015 OpenFOAM Foundation
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
    Foam::myIncompressibleMyTurbulenceModel

Description
    Abstract base class for turbulence models (RAS, LES and laminar).

SourceFiles
    myIncompressibleMyTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef myIncompressibleMyTurbulenceModel_H
#define myIncompressibleMyTurbulenceModel_H

#include "myTurbulenceModel.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                Class myIncompressibleMyTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class myIncompressibleMyTurbulenceModel
:
    public myTurbulenceModel
{

protected:

    // Protected data

        geometricOneField rho_;


    // Protected member functions

        //- ***HGW Temporary function to be removed when the run-time selectable
        //  thermal transport layer is complete
        virtual void correctNut()
        {}


private:

    // Private Member Functions

        //- No copy construct
        myIncompressibleMyTurbulenceModel
        (
            const myIncompressibleMyTurbulenceModel&
        ) = delete;

        //- No copy assignment
        void operator=(const myIncompressibleMyTurbulenceModel&) = delete;


public:

    //- Runtime type information
    TypeName("myIncompressibleMyTurbulenceModel");


    // Constructors

        //- Construct from components
        myIncompressibleMyTurbulenceModel
        (
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const word& propertiesName
        );


    //- Destructor
    virtual ~myIncompressibleMyTurbulenceModel()
    {}


    // Member Functions

        //- Return the laminar dynamic viscosity
        virtual tmp<volScalarField> mu() const;

        //- Return the laminar dynamic viscosity on patch
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Return the turbulence dynamic viscosity
        virtual tmp<volScalarField> mut() const;

        //- Return the turbulence dynamic viscosity on patch
        virtual tmp<scalarField> mut(const label patchi) const;

        //- Return the effective dynamic viscosity
        virtual tmp<volScalarField> muEff() const;

        //- Return the effective dynamic viscosity on patch
        virtual tmp<scalarField> muEff(const label patchi) const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const = 0;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
