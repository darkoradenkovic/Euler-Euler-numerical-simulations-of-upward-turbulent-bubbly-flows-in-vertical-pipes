/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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

Class
    Foam::myOmegaWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall constraint on turbulence specific
    dissipation, omega for both low and high Reynolds number turbulence models.

    The near-wall omega may be either blended between the viscous region and
    logarithmic region values using:

        \f[
            \omega = sqrt(\omega_{vis}^2 + \omega_{log}^2)
        \f]

    where

    \vartable
        \omega_{vis} | omega in viscous region
        \omega_{log} | omega in logarithmic region
    \endvartable

    see eq.(15) of:
    \verbatim
        Menter, F., Esch, T.
        "Elements of Industrial Heat Transfer Prediction"
        16th Brazilian Congress of Mechanical Engineering (COBEM),
        Nov. 2001
    \endverbatim

    or switched between these values based on the laminar-to-turbulent y+ value
    derived from kappa and E.  Recent tests have shown that the standard
    switching method provides more accurate results for 10 < y+ < 30 when used
    with high Reynolds number wall-functions and both methods provide accurate
    results when used with continuous wall-functions.  Based on this the
    standard switching method is used by default.

Usage
    \table
        Property     | Description             | Required    | Default value
        Cmu          | Model coefficient       | no          | 0.09
        kappa        | von Karman constant     | no          | 0.41
        E            | Model coefficient       | no          | 9.8
        beta1        | Model coefficient       | no          | 0.075
        blended      | Blending switch         | no          | false
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            myOmegaWallFunction;
    }
    \endverbatim

See also
    Foam::fixedInternalValueFvPatchField
    Foam::myEpsilonWallFunctionFvPatchScalarField

SourceFiles
    myOmegaWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef myOmegaWallFunctionFvPatchScalarField_H
#define myOmegaWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myTurbulenceModel;

/*---------------------------------------------------------------------------*\
             Class myOmegaWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class myOmegaWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

        //- Tolerance used in weighted calculations
        static scalar tolerance_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;

        //- beta1 coefficient
        scalar beta1_;

        //- beta1 coefficient
        Switch blended_;

        //- y+ at the edge of the laminar sublayer
        scalar yPlusLam_;

        //- Local copy of turbulence G field
        scalarField G_;

        //- Local copy of turbulence omega field
        scalarField omega_;

        //- Initialised flag
        bool initialised_;

        //- Master patch ID
        label master_;

        //- List of averaging corner weights
        List<List<scalar>> cornerWeights_;


    // Protected Member Functions

        //- Check the type of the patch
        virtual void checkType();

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;

        //- Set the master patch - master is responsible for updating all
        //  wall function patches
        virtual void setMaster();

        //- Create the averaging weights for cells which are bounded by
        //  multiple wall function faces
        virtual void createAveragingWeights();

        //- Helper function to return non-const access to an omega patch
        virtual myOmegaWallFunctionFvPatchScalarField& omegaPatch
        (
            const label patchi
        );

        //- Main driver to calculate the turbulence fields
        virtual void calculateTurbulenceFields
        (
            const myTurbulenceModel& turbulence,
            scalarField& G0,
            scalarField& omega0
        );

        //- Calculate the omega and G
        virtual void calculate
        (
            const myTurbulenceModel& turbulence,
            const List<scalar>& cornerWeights,
            const fvPatch& patch,
            scalarField& G,
            scalarField& omega
        );

        //- Return non-const access to the master patch ID
        virtual label& master()
        {
            return master_;
        }


public:

    //- Runtime type information
    TypeName("myOmegaWallFunction");


    // Constructors

        //- Construct from patch and internal field
        myOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myOmegaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // myOmegaWallFunctionFvPatchScalarField
        //  onto a new patch
        myOmegaWallFunctionFvPatchScalarField
        (
            const myOmegaWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myOmegaWallFunctionFvPatchScalarField
        (
            const myOmegaWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new myOmegaWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myOmegaWallFunctionFvPatchScalarField
        (
            const myOmegaWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new myOmegaWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return non-const access to the master's G field
            scalarField& G(bool init = false);

            //- Return non-const access to the master's omega field
            scalarField& omega(bool init = false);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Update the coefficients associated with the patch field
            virtual void updateWeightedCoeffs(const scalarField& weights);

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);

            //- Manipulate matrix with given weights
            virtual void manipulateMatrix
            (
                fvMatrix<scalar>& matrix,
                const scalarField& weights
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
