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
    Foam::phaseCompressibleMyTurbulenceModel

Description
    Forward declaration of typedef for phaseCompressibleMyTurbulenceModel

\*---------------------------------------------------------------------------*/

#ifndef phaseCompressibleMyTurbulenceModelFwd_H
#define phaseCompressibleMyTurbulenceModelFwd_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class myPhaseModel;

    template<class TransportModel>
    class PhaseCompressibleMyTurbulenceModel;

    template<class BasicMyTurbulenceModel>
    class ThermalDiffusivity;

      typedef ThermalDiffusivity<PhaseCompressibleMyTurbulenceModel<myPhaseModel>>
          phaseCompressibleMyTurbulenceModel;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
