/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2011 OpenCFD Ltd.
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

#include "turbulentMixingLengthDissipationRateInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "myTurbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    mixingLength_(0.0),
    kName_("k")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    mixingLength_(ptf.mixingLength_),
    kName_(ptf.kName_)
{}


turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    mixingLength_(dict.get<scalar>("mixingLength")),
    kName_(dict.lookupOrDefault<word>("k", "k"))
{
    this->phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    mixingLength_(ptf.mixingLength_),
    kName_(ptf.kName_)
{}


turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    mixingLength_(ptf.mixingLength_),
    kName_(ptf.kName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentMixingLengthDissipationRateInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup Cmu corresponding to the turbulence model selected
    const myTurbulenceModel& turbModel = db().lookupObject<myTurbulenceModel>
    (
        IOobject::groupName
        (
            myTurbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalar Cmu =
        turbModel.coeffDict().lookupOrDefault<scalar>("Cmu", 0.09);

    const scalar Cmu75 = pow(Cmu, 0.75);

    const fvPatchScalarField& kp =
        patch().lookupPatchField<volScalarField, scalar>(kName_);

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    this->refValue() = Cmu75*kp*sqrt(kp)/mixingLength_;
    this->valueFraction() = 1.0 - pos0(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void turbulentMixingLengthDissipationRateInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("mixingLength", mixingLength_);
    os.writeEntry("phi", this->phiName_);
    os.writeEntry("k", kName_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentMixingLengthDissipationRateInletFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
