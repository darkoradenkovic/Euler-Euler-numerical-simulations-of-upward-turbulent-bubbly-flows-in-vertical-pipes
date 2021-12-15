/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2015-2016 OpenFOAM Foundation
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

Typedef
    Foam::compressible::myTurbulenceModel

Typedef
    Foam::compressible::RASModel

Typedef
    Foam::compressible::LESModel

Description
    Typedefs for turbulence, RAS and LES models for compressible flow
    based on the standard laminar transport package.

SourceFiles
    turbulentMyFluidThermoModel.C
    turbulentMyFluidThermoModels.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentMyFluidThermoModel_H
#define turbulentMyFluidThermoModel_H

#include "CompressibleMyTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "EddyDiffusivity.H"
#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "myFluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace compressible
    {
        typedef ThermalDiffusivity<CompressibleMyTurbulenceModel<myFluidThermo>>
            myTurbulenceModel;

        typedef laminarModel<myTurbulenceModel> laminarModel;
        typedef RASModel<EddyDiffusivity<myTurbulenceModel>> RASModel;
        typedef LESModel<EddyDiffusivity<myTurbulenceModel>> LESModel;

        template<class BasicCompressibleMyTurbulenceModel>
        autoPtr<BasicCompressibleMyTurbulenceModel> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const typename BasicCompressibleMyTurbulenceModel::transportModel&
                transport,
            const word& propertiesName = myTurbulenceModel::propertiesName
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentMyFluidThermoModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
