/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "PhaseCompressibleMyTurbulenceModel.H"
#include "myPhaseModel.H"
#include "myTwoPhaseSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "makeMyTurbulenceModel.H"

#include "ThermalDiffusivity.H"
#include "EddyDiffusivity.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeMyTurbulenceModelTypes
(
    volScalarField,
    volScalarField,
    myCompressibleMyTurbulenceModel,
    PhaseCompressibleMyTurbulenceModel,
    ThermalDiffusivity,
    myPhaseModel
);

makeBaseMyTurbulenceModel
(
    volScalarField,
    volScalarField,
    myCompressibleMyTurbulenceModel,
    PhaseCompressibleMyTurbulenceModel,
    ThermalDiffusivity,
    myPhaseModel
);

#define makeLaminarModel(Type)                                                 \
    makeTemplatedLaminarModel                                                  \
    (myPhaseModelPhaseCompressibleMyTurbulenceModel, laminar, Type)

#define makeRASModel(Type)                                                     \
    makeTemplatedMyTurbulenceModel                                               \
    (myPhaseModelPhaseCompressibleMyTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedMyTurbulenceModel                                               \
    (myPhaseModelPhaseCompressibleMyTurbulenceModel, LES, Type)

#include "Stokes.H"
makeLaminarModel(Stokes);

#include "kEpsilon.H"
makeRASModel(kEpsilon);

#include "kOmegaSST.H"
makeRASModel(kOmegaSST);

#include "kOmegaSSTSato.H"
makeRASModel(kOmegaSSTSato);

#include "SSG.H"
makeRASModel (SSG);

#include "LRR.H"
makeRASModel (LRR);

#include "LaunderSharmaKE.H"
makeRASModel (LaunderSharmaKE);

#include "LaunderSharmaKECF.H"
makeRASModel (LaunderSharmaKECF);

#include "mixtureKEpsilon.H"
makeRASModel(mixtureKEpsilon);

#include "LaheyKEpsilon.H"
makeRASModel(LaheyKEpsilon);

#include "continuousGasKEpsilon.H"
makeRASModel(continuousGasKEpsilon);

#include "EBRSM.H"
makeRASModel(EBRSM);

#include "Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "kEqn.H"
makeLESModel(kEqn);

#include "SmagorinskyZhang.H"
makeLESModel(SmagorinskyZhang);

#include "NicenoKEqn.H"
makeLESModel(NicenoKEqn);

#include "continuousGasKEqn.H"
makeLESModel(continuousGasKEqn);

#include "kineticTheoryModel.H"
makeMyTurbulenceModel
(myPhaseModelPhaseCompressibleMyTurbulenceModel, RAS, kineticTheoryModel);

#include "phasePressureModel.H"
makeMyTurbulenceModel
(myPhaseModelPhaseCompressibleMyTurbulenceModel, RAS, phasePressureModel);

// ************************************************************************* //
