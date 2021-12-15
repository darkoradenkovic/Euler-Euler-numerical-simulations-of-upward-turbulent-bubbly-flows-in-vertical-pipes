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

Class
    Foam::myBlendedInterfacialModel

Description

SourceFiles
    myBlendedInterfacialModel.C

\*---------------------------------------------------------------------------*/

#ifndef myBlendedInterfacialModel_H
#define myBlendedInterfacialModel_H

#include "blendingMethod.H"
#include "myPhasePair.H"
#include "orderedPhasePair.H"

#include "geometricZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class myBlendedInterfacialModel Declaration
\*---------------------------------------------------------------------------*/

template<class modelType>
class myBlendedInterfacialModel
{
    // Private data

        //- Unordered phase pair
        const myPhasePair& pair_;

        //- Ordered phase pair for dispersed phase 1 in continuous phase 2
        const orderedPhasePair& pair1In2_;

        //- Ordered phase pair for dispersed phase 2 in continuous phase 1
        const orderedPhasePair& pair2In1_;

        //- Model for region with no obvious dispersed phase
        autoPtr<modelType> model_;

        //- Model for dispersed phase 1 in continuous phase 2
        autoPtr<modelType> model1In2_;

        //- Model for dispersed phase 2 in continuous phase 1
        autoPtr<modelType> model2In1_;

        //- Blending model
        const blendingMethod& blending_;

        //- If true set coefficients and forces to 0 at fixed-flux BCs
        bool correctFixedFluxBCs_;


    // Private Member Functions

        //- No copy construct
        myBlendedInterfacialModel
        (
            const myBlendedInterfacialModel&
        ) = delete;

        //- No copy assignment
        void operator=(const myBlendedInterfacialModel<modelType>&) = delete;

        //- Correct coeff/value on fixed flux boundary conditions
        template<class GeometricField>
        void correctFixedFluxBCs(GeometricField& field) const;


public:

    // Constructors

        //- Construct from the model table, dictionary and pairs
        myBlendedInterfacialModel
        (
            const myPhasePair::dictTable& modelTable,
            const blendingMethod& blending,
            const myPhasePair& pair,
            const orderedPhasePair& pair1In2,
            const orderedPhasePair& pair2In1,
            const bool correctFixedFluxBCs = true
        );


    //- Destructor
    ~myBlendedInterfacialModel();


    // Member Functions

        //- Return true if a model is specified for the supplied phase
        bool hasModel(const myPhaseModel& phase) const;

        //- Return the model for the supplied phase
        const modelType& myPhaseModel(const myPhaseModel& phase) const;

        //- Return the blended force coefficient
        tmp<volScalarField> K() const;

        //- Return the face blended force coefficient
        tmp<surfaceScalarField> Kf() const;

        //- Return the blended force
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> F() const;
        
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> FLam() const;
        
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> FTurb() const;

        //- Return the face blended force
        tmp<surfaceScalarField> Ff() const;

        //- Return the blended diffusivity
        tmp<volScalarField> D() const;
        
        //- Return the blended diffusivity
        tmp<volScalarField> DLam() const;
        
              //- Return the blended diffusivity
        tmp<volScalarField> DTurb() const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "myBlendedInterfacialModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
