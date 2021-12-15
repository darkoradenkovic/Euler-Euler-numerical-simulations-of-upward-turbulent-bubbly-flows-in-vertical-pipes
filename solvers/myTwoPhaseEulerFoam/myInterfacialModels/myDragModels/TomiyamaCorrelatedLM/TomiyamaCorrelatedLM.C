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

#include "TomiyamaCorrelatedLM.H"
#include "myPhasePair.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myDragModels
{
    defineTypeNameAndDebug(TomiyamaCorrelatedLM, 0);
    addToRunTimeSelectionTable(myDragModel, TomiyamaCorrelatedLM, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myDragModels::TomiyamaCorrelatedLM::TomiyamaCorrelatedLM
(
    const dictionary& dict,
    const myPhasePair& pair,
    const bool registerObject
)
:
    myDragModel(dict, pair, registerObject),
    A_("A", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myDragModels::TomiyamaCorrelatedLM::~TomiyamaCorrelatedLM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::myDragModels::TomiyamaCorrelatedLM::CdRe() const
{
    volScalarField Re(pair_.Re());
    volScalarField Eo(pair_.Eo());

    volScalarField Sr
    (
        sqr(pair_.dispersed().d())
       /(
            Re
           *pair_.continuous().nu()
        )
       *mag(fvc::grad(pair_.continuous().U()))
    );
    
    return
        (max
        (
            A_
           *min
            (
                (1 + 0.15*pow(Re, 0.687)),
                scalar(3)
            ),
            8*Eo*Re/(3*Eo + 12))*(1.0+0.55*sqr(Sr))
        );

}


// ************************************************************************* //
