/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "SyamlalOBrien.H"
#include "myPhasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myDragModels
{
    defineTypeNameAndDebug(SyamlalOBrien, 0);
    addToRunTimeSelectionTable(myDragModel, SyamlalOBrien, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myDragModels::SyamlalOBrien::SyamlalOBrien
(
    const dictionary& dict,
    const myPhasePair& pair,
    const bool registerObject
)
:
    myDragModel(dict, pair, registerObject)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myDragModels::SyamlalOBrien::~SyamlalOBrien()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::myDragModels::SyamlalOBrien::CdRe() const
{
    volScalarField alpha2
    (
        max(scalar(1) - pair_.dispersed(), pair_.continuous().residualAlpha())
    );

    volScalarField A(pow(alpha2, 4.14));
    volScalarField B
    (
        neg(alpha2 - 0.85)*(0.8*pow(alpha2, 1.28))
      + pos0(alpha2 - 0.85)*(pow(alpha2, 2.65))
    );
    volScalarField Re(pair_.Re());
    volScalarField Vr
    (
        0.5
       *(
            A - 0.06*Re + sqrt(sqr(0.06*Re) + 0.12*Re*(2.0*B - A) + sqr(A))
        )
    );
    volScalarField CdsRe(sqr(0.63*sqrt(Re) + 4.8*sqrt(Vr)));

    return
        CdsRe
       *max(pair_.continuous(), pair_.continuous().residualAlpha())
       /sqr(Vr);
}


// ************************************************************************* //
