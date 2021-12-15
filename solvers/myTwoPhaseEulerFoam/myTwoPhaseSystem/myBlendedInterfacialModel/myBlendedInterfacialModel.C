/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "myBlendedInterfacialModel.H"
#include "fixedValueFvsPatchFields.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class modelType>
template<class GeometricField>
void Foam::myBlendedInterfacialModel<modelType>::correctFixedFluxBCs
(
    GeometricField& field
) const
{
    typename GeometricField::Boundary& fieldBf =
        field.boundaryFieldRef();

    forAll(pair_.phase1().phi().boundaryField(), patchi)
    {
        if
        (
            isA<fixedValueFvsPatchScalarField>
            (
                pair_.phase1().phi().boundaryField()[patchi]
            )
        )
        {
            fieldBf[patchi] = Zero;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class modelType>
Foam::myBlendedInterfacialModel<modelType>::myBlendedInterfacialModel
(
    const myPhasePair::dictTable& modelTable,
    const blendingMethod& blending,
    const myPhasePair& pair,
    const orderedPhasePair& pair1In2,
    const orderedPhasePair& pair2In1,
    const bool correctFixedFluxBCs
)
:
    pair_(pair),
    pair1In2_(pair1In2),
    pair2In1_(pair2In1),
    blending_(blending),
    correctFixedFluxBCs_(correctFixedFluxBCs)
{
    if (modelTable.found(pair_))
    {
        model_.reset
        (
            modelType::New
            (
                modelTable[pair_],
                pair_
            ).ptr()
        );
    }

    if (modelTable.found(pair1In2_))
    {
        model1In2_.reset
        (
            modelType::New
            (
                modelTable[pair1In2_],
                pair1In2_
            ).ptr()
        );
    }

    if (modelTable.found(pair2In1_))
    {
        model2In1_.reset
        (
            modelType::New
            (
                modelTable[pair2In1_],
                pair2In1_
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class modelType>
Foam::myBlendedInterfacialModel<modelType>::~myBlendedInterfacialModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class modelType>
Foam::tmp<Foam::volScalarField>
Foam::myBlendedInterfacialModel<modelType>::K() const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    tmp<volScalarField> x
    (
        new volScalarField
        (
            IOobject
            (
                modelType::typeName + ":K",
                pair_.phase1().mesh().time().timeName(),
                pair_.phase1().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pair_.phase1().mesh(),
            dimensionedScalar(modelType::dimK, Zero)
        )
    );

    if (model_.valid())
    {
        x.ref() += model_->K()*(f1() - f2());
    }

    if (model1In2_.valid())
    {
        x.ref() += model1In2_->K()*(1 - f1);
    }

    if (model2In1_.valid())
    {
        x.ref() += model2In1_->K()*f2;
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class modelType>
Foam::tmp<Foam::surfaceScalarField>
Foam::myBlendedInterfacialModel<modelType>::Kf() const
{
    tmp<surfaceScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = fvc::interpolate
        (
            blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed())
        );
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = fvc::interpolate
        (
            blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed())
        );
    }

    tmp<surfaceScalarField> x
    (
        new surfaceScalarField
        (
            IOobject
            (
                modelType::typeName + ":Kf",
                pair_.phase1().mesh().time().timeName(),
                pair_.phase1().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pair_.phase1().mesh(),
            dimensionedScalar(modelType::dimK, Zero)
        )
    );

    if (model_.valid())
    {
        x.ref() += model_->Kf()*(f1() - f2());
    }

    if (model1In2_.valid())
    {
        x.ref() += model1In2_->Kf()*(1 - f1);
    }

    if (model2In1_.valid())
    {
        x.ref() += model2In1_->Kf()*f2;
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class modelType>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::myBlendedInterfacialModel<modelType>::F() const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    auto x = tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            modelType::typeName + ":F",
            pair_.phase1().mesh().time().timeName(),
            pair_.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pair_.phase1().mesh(),
        dimensioned<Type>(modelType::dimF, Zero)
    );

    if (model_.valid())
    {
        x.ref() += model_->F()*(f1() - f2());
    }

    if (model1In2_.valid())
    {
        x.ref() += model1In2_->F()*(1 - f1);
    }

    if (model2In1_.valid())
    {
        x.ref() -= model2In1_->F()*f2; // note : subtraction
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class modelType>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::myBlendedInterfacialModel<modelType>::FLam() const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    auto x = tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            modelType::typeName + ":FLam",
            pair_.phase1().mesh().time().timeName(),
            pair_.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pair_.phase1().mesh(),
        dimensioned<Type>(modelType::dimF, Zero)
    );

    if (model_.valid())
    {
        x.ref() += model_->FLam()*(f1() - f2());
    }

    if (model1In2_.valid())
    {
        x.ref() += model1In2_->FLam()*(1 - f1);
    }

    if (model2In1_.valid())
    {
        x.ref() -= model2In1_->FLam()*f2; // note : subtraction
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class modelType>
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::myBlendedInterfacialModel<modelType>::FTurb() const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    auto x = tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        IOobject
        (
            modelType::typeName + ":FTurb",
            pair_.phase1().mesh().time().timeName(),
            pair_.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pair_.phase1().mesh(),
        dimensioned<Type>(modelType::dimF, Zero)
    );

    if (model_.valid())
    {
        x.ref() += model_->FTurb()*(f1() - f2());
    }

    if (model1In2_.valid())
    {
        x.ref() += model1In2_->FTurb()*(1 - f1);
    }

    if (model2In1_.valid())
    {
        x.ref() -= model2In1_->FTurb()*f2; // note : subtraction
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}

template<class modelType>
Foam::tmp<Foam::surfaceScalarField>
Foam::myBlendedInterfacialModel<modelType>::Ff() const
{
    tmp<surfaceScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = fvc::interpolate
        (
            blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed())
        );
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = fvc::interpolate
        (
            blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed())
        );
    }

    auto x = tmp<surfaceScalarField>::New
    (
        IOobject
        (
            modelType::typeName + ":Ff",
            pair_.phase1().mesh().time().timeName(),
            pair_.phase1().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pair_.phase1().mesh(),
        dimensionedScalar(modelType::dimF*dimArea, Zero)
    );

    x.ref().setOriented();

    if (model_.valid())
    {
        x.ref() += model_->Ff()*(f1() - f2());
    }

    if (model1In2_.valid())
    {
        x.ref() += model1In2_->Ff()*(1 - f1);
    }

    if (model2In1_.valid())
    {
        x.ref() -= model2In1_->Ff()*f2; // note : subtraction
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class modelType>
Foam::tmp<Foam::volScalarField>
Foam::myBlendedInterfacialModel<modelType>::D() const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    tmp<volScalarField> x
    (
        new volScalarField
        (
            IOobject
            (
                modelType::typeName + ":D",
                pair_.phase1().mesh().time().timeName(),
                pair_.phase1().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pair_.phase1().mesh(),
            dimensionedScalar(modelType::dimD, Zero)
        )
    );

    if (model_.valid())
    {
        x.ref() += model_->D()*(f1() - f2());
    }

    if (model1In2_.valid())
    {
        x.ref() += model1In2_->D()*(1 - f1);
    }

    if (model2In1_.valid())
    {
        x.ref() += model2In1_->D()*f2;
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}

template<class modelType>
Foam::tmp<Foam::volScalarField>
Foam::myBlendedInterfacialModel<modelType>::DLam() const
{
    tmp<volScalarField> f1, f2;

    if (model_.valid() || model1In2_.valid())
    {
        f1 = blending_.f1(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    if (model_.valid() || model2In1_.valid())
    {
        f2 = blending_.f2(pair1In2_.dispersed(), pair2In1_.dispersed());
    }

    tmp<volScalarField> x
    (
        new volScalarField
        (
            IOobject
            (
                modelType::typeName + ":DLam",
                pair_.phase1().mesh().time().timeName(),
                pair_.phase1().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pair_.phase1().mesh(),
            dimensionedScalar(modelType::dimD, Zero)
        )
    );

    if (model_.valid())
    {
        x.ref() += model_->DLam()*(f1() - f2());
    }

    if (model1In2_.valid())
    {
        x.ref() += model1In2_->DLam()*(1 - f1);
    }

    if (model2In1_.valid())
    {
        x.ref() += model2In1_->DLam()*f2;
    }

    if
    (
        correctFixedFluxBCs_
     && (model_.valid() || model1In2_.valid() || model2In1_.valid())
    )
    {
        correctFixedFluxBCs(x.ref());
    }

    return x;
}


template<class modelType>
bool Foam::myBlendedInterfacialModel<modelType>::hasModel
(
    const class myPhaseModel& phase
) const
{
    return
    (
        &phase == &(pair_.phase1())
      ? model1In2_.valid()
      : model2In1_.valid()
    );
}


template<class modelType>
const modelType& Foam::myBlendedInterfacialModel<modelType>::myPhaseModel
(
    const class myPhaseModel& phase
) const
{
    return &phase == &(pair_.phase1()) ? *model1In2_ : *model2In1_;
}


// ************************************************************************* //
