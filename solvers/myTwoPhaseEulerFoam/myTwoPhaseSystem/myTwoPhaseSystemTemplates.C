/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "myBlendedInterfacialModel.H"
#include "myDragModel.H"
#include "myVirtualMassModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class modelType>
const modelType& myTwoPhaseSystem::lookupSubModel
(
    const myPhasePair& key
) const
{
    return
        mesh().lookupObject<modelType>
        (
            IOobject::groupName(modelType::typeName, key.name())
        );
}


template<>
inline const myDragModel& myTwoPhaseSystem::lookupSubModel<myDragModel>
(
    const myPhaseModel& dispersed,
    const myPhaseModel& continuous
) const
{
    return drag_->myPhaseModel(dispersed);
}


template<>
inline const myVirtualMassModel& myTwoPhaseSystem::lookupSubModel<myVirtualMassModel>
(
    const myPhaseModel& dispersed,
    const myPhaseModel& continuous
) const
{
    return virtualMass_->myPhaseModel(dispersed);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
