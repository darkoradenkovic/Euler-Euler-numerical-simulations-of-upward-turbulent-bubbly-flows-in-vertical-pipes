/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

Class
    Foam::myKqRWallFunctionFvPatchField

Group
    grpWallFunctions

Description
    This boundary condition provides a suitable condition for turbulence
    \c k, \c q, and \c R fields for the case of high Reynolds number flow using
    wall functions.

    It is a simple wrapper around the zero-gradient condition.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            myKqRWallFunction;
    }
    \endverbatim

See also
    Foam::zeroGradientFvPatchField

SourceFiles
    myKqRWallFunctionFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef myKqRWallFunctionFvPatchField_H
#define myKqRWallFunctionFvPatchField_H

#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class myKqRWallFunctionFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class myKqRWallFunctionFvPatchField
:
    public zeroGradientFvPatchField<Type>
{

    // Private Member Functions

        //- Check the type of the patch
        void checkType();


public:

    //- Runtime type information
    TypeName("myKqRWallFunction");


    // Constructors

        //- Construct from patch and internal field
        myKqRWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        myKqRWallFunctionFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // myKqRWallFunctionFvPatchField
        //  onto a new patch
        myKqRWallFunctionFvPatchField
        (
            const myKqRWallFunctionFvPatchField&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        myKqRWallFunctionFvPatchField
        (
            const myKqRWallFunctionFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new myKqRWallFunctionFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        myKqRWallFunctionFvPatchField
        (
            const myKqRWallFunctionFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new myKqRWallFunctionFvPatchField(*this, iF)
            );
        }


    // Member functions

        //- Evaluate the patchField
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "myKqRWallFunctionFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
