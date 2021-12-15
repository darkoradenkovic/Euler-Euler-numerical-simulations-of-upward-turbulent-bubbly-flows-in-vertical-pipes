/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "myNutKWallFunctionFvPatchScalarField.H"
#include "myTurbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> myNutKWallFunctionFvPatchScalarField::calcNut() const
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
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(patch().size(), Zero));
    scalarField& nutw = tnutw.ref();

    forAll(nutw, facei)
    {
        label celli = patch().faceCells()[facei];

        scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];

        if (yPlus > yPlusLam_)
        {
            nutw[facei] = nuw[facei]*(yPlus*kappa_/log(E_*yPlus) - 1.0);
        }
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myNutKWallFunctionFvPatchScalarField::myNutKWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    myNutWallFunctionFvPatchScalarField(p, iF)
{}


myNutKWallFunctionFvPatchScalarField::myNutKWallFunctionFvPatchScalarField
(
    const myNutKWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    myNutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


myNutKWallFunctionFvPatchScalarField::myNutKWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    myNutWallFunctionFvPatchScalarField(p, iF, dict)
{}


myNutKWallFunctionFvPatchScalarField::myNutKWallFunctionFvPatchScalarField
(
    const myNutKWallFunctionFvPatchScalarField& wfpsf
)
:
    myNutWallFunctionFvPatchScalarField(wfpsf)
{}


myNutKWallFunctionFvPatchScalarField::myNutKWallFunctionFvPatchScalarField
(
    const myNutKWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    myNutWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> myNutKWallFunctionFvPatchScalarField::yPlus() const
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

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    tmp<scalarField> kwc = k.boundaryField()[patchi].patchInternalField();
    tmp<scalarField> nuw = turbModel.nu(patchi);

    return pow025(Cmu_)*y*sqrt(kwc)/nuw;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    myNutKWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
