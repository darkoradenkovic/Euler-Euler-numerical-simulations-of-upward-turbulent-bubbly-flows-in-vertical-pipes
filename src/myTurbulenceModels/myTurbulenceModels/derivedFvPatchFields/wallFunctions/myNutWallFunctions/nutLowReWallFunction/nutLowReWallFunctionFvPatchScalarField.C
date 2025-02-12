/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "nutLowReWallFunctionFvPatchScalarField.H"
#include "myTurbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutLowReWallFunctionFvPatchScalarField::calcNut() const
{
    return tmp<scalarField>::New(patch().size(), Zero);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutLowReWallFunctionFvPatchScalarField::nutLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    myNutWallFunctionFvPatchScalarField(p, iF)
{}


nutLowReWallFunctionFvPatchScalarField::nutLowReWallFunctionFvPatchScalarField
(
    const nutLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    myNutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


nutLowReWallFunctionFvPatchScalarField::nutLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    myNutWallFunctionFvPatchScalarField(p, iF, dict)
{}


nutLowReWallFunctionFvPatchScalarField::nutLowReWallFunctionFvPatchScalarField
(
    const nutLowReWallFunctionFvPatchScalarField& nlrwfpsf
)
:
    myNutWallFunctionFvPatchScalarField(nlrwfpsf)
{}


nutLowReWallFunctionFvPatchScalarField::nutLowReWallFunctionFvPatchScalarField
(
    const nutLowReWallFunctionFvPatchScalarField& nlrwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    myNutWallFunctionFvPatchScalarField(nlrwfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutLowReWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();
    const myTurbulenceModel& turbModel = db().lookupObject<myTurbulenceModel>
    (
        IOobject::groupName
        (
            myTurbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];

    return y*sqrt(nuw*mag(Uw.snGrad()))/nuw;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutLowReWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
