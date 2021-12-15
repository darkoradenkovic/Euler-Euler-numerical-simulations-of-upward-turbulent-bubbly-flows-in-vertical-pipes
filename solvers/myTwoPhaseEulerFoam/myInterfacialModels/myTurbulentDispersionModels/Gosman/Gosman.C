/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2015 OpenFOAM Foundation
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

#include "Gosman.H"
#include "myPhasePair.H"
#include "PhaseCompressibleMyTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

#include "myDragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myTurbulentDispersionModels
{
    defineTypeNameAndDebug(Gosman, 0);
    addToRunTimeSelectionTable
    (
        myTurbulentDispersionModel,
        Gosman,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myTurbulentDispersionModels::Gosman::Gosman
(
    const dictionary& dict,
    const myPhasePair& pair
)
:
    myTurbulentDispersionModel(dict, pair),
    sigma_("sigma", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myTurbulentDispersionModels::Gosman::~Gosman()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::myTurbulentDispersionModels::Gosman::D() const
{
    const fvMesh& mesh(pair_.phase1().mesh());
    const myDragModel&
        drag
        (
            mesh.lookupObject<myDragModel>
            (
                IOobject::groupName(myDragModel::typeName, pair_.name())
            )
        );

    return
        0.75
       *drag.CdRe()
       *pair_.dispersed()
       *pair_.continuous().nu()
       *pair_.continuous().turbulence().nut()
       /(
            sigma_
           *sqr(pair_.dispersed().d())
        )
       *pair_.continuous().rho();
}


// ************************************************************************* //
