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

#include "completeDispersion.H"
#include "myPhasePair.H"
#include "PhaseCompressibleMyTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

#include "myDragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myTurbulentDispersionModels
{
    defineTypeNameAndDebug(completeDispersion, 0);
    addToRunTimeSelectionTable
    (
        myTurbulentDispersionModel,
        completeDispersion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myTurbulentDispersionModels::completeDispersion::completeDispersion
(
    const dictionary& dict,
    const myPhasePair& pair
)
:
    myTurbulentDispersionModel(dict, pair),
    c_("C", dimless, dict),
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

Foam::myTurbulentDispersionModels::completeDispersion::~completeDispersion()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::myTurbulentDispersionModels::completeDispersion::D() const
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
        
  const volVectorField& DragForce =
            mesh.lookupObject<volVectorField>("DragForce");
        
        
    return
       0.33333*pair_.dispersed()*mag(DragForce)*pair_.dispersed().d()*(2*c_-1.0)*(1.0-4.0*pair_.dispersed()) +
       0.75
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
        );

}





Foam::tmp<Foam::volScalarField>
Foam::myTurbulentDispersionModels::completeDispersion::DLam() const
{
    const fvMesh& mesh(pair_.phase1().mesh());
    
     
  const volVectorField& DragForce =
            mesh.lookupObject<volVectorField>("DragForce");
        
        
    return
    0.33333*pair_.dispersed()*mag(DragForce)*pair_.dispersed().d()*(2*c_-1.0)*(1.0-4.0*pair_.dispersed());

}


Foam::tmp<Foam::volScalarField>
Foam::myTurbulentDispersionModels::completeDispersion::DTurb() const
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
        );

}


// ************************************************************************* //
