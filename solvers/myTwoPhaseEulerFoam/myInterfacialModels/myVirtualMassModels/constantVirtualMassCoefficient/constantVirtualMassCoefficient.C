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

#include "constantVirtualMassCoefficient.H"
#include "myPhasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace myVirtualMassModels
{
    defineTypeNameAndDebug(constantVirtualMassCoefficient, 0);
    addToRunTimeSelectionTable
    (
        myVirtualMassModel,
        constantVirtualMassCoefficient,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myVirtualMassModels::constantVirtualMassCoefficient::
constantVirtualMassCoefficient
(
    const dictionary& dict,
    const myPhasePair& pair,
    const bool registerObject
)
:
    myVirtualMassModel(dict, pair, registerObject),
    Cvm_("Cvm", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::myVirtualMassModels::constantVirtualMassCoefficient::
~constantVirtualMassCoefficient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::myVirtualMassModels::constantVirtualMassCoefficient::Cvm() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Cvm",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            Cvm_
        )
    );
}


// ************************************************************************* //
