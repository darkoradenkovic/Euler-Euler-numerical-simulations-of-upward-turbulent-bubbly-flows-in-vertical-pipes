/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2017 OpenFOAM Foundation
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

#include "Lubchenko.H"
#include "myPhasePair.H"
#include "PhaseCompressibleMyTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

#include "myDragModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myWallLubricationModels
{
    defineTypeNameAndDebug(Lubchenko, 0);
    addToRunTimeSelectionTable
    (
        myWallLubricationModel,
        Lubchenko,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myWallLubricationModels::Lubchenko::Lubchenko
(
    const dictionary& dict,
    const myPhasePair& pair
)
:
    myWallLubricationModel(dict, pair),
    sigma_("sigma", dimless, dict),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        pair_.dispersed().residualAlpha().value(),
        dict
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myWallLubricationModels::Lubchenko::~Lubchenko()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::myWallLubricationModels::Lubchenko::Fi() const
{
    const volVectorField& n(nWall());
    const volScalarField& y(yWall());  
        
    
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
        (
        neg(y/pair_.dispersed().d()-0.5)*1.0  
       +pos0(y/pair_.dispersed().d()-0.5)*0 
        )
        *(0.75)
       *drag.CdRe()
       *pair_.continuous().nu()
       *pair_.continuous().turbulence().nut()
       /(
            sigma_
           *sqr(pair_.dispersed().d())
        )
       *pair_.continuous().rho()
       *pair_.dispersed()
       *(
           1.0/max(pair_.dispersed(), residualAlpha_)
         + 1.0/max(pair_.continuous(), residualAlpha_)
        )
       *1.0/y*(pair_.dispersed().d()-2.0*y)/(pair_.dispersed().d()-y)*n;
}


// ************************************************************************* //
