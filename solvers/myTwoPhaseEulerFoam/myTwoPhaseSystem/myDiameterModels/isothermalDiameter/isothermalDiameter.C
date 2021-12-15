/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "isothermalDiameter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myDiameterModels
{
    defineTypeNameAndDebug(isothermal, 0);

    addToRunTimeSelectionTable
    (
        myDiameterModel,
        isothermal,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myDiameterModels::isothermal::isothermal
(
    const dictionary& diameterProperties,
    const myPhaseModel& phase
)
:
    myDiameterModel(diameterProperties, phase),
    d0_("d0", dimLength, diameterProperties_),
    p0_("p0", dimPressure, diameterProperties_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myDiameterModels::isothermal::~isothermal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::myDiameterModels::isothermal::d() const
{
    const volScalarField& p = phase_.U().db().lookupObject<volScalarField>
    (
        "p"
    );

    return d0_*cbrt(p0_/p);
}


bool Foam::myDiameterModels::isothermal::read(const dictionary& phaseProperties)
{
    myDiameterModel::read(phaseProperties);

    diameterProperties_.readEntry("d0", d0_);
    diameterProperties_.readEntry("p0", p0_);

    return true;
}


// ************************************************************************* //
