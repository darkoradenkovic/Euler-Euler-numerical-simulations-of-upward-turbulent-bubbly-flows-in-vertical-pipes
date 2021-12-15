/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::myPhasePair

Description

SourceFiles
    myPhasePair.C

\*---------------------------------------------------------------------------*/

#ifndef myPhasePair_H
#define myPhasePair_H

#include "myPhaseModel.H"
#include "myPhasePairKey.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class myPhasePair Declaration
\*---------------------------------------------------------------------------*/

class myPhasePair
:
    public myPhasePairKey
{
public:

    // Hash table types

        //- Dictionary hash table
        typedef HashTable<dictionary, myPhasePairKey, myPhasePairKey::hash>
            dictTable;

        //- Scalar hash table
        typedef HashTable<scalar, myPhasePairKey, myPhasePairKey::hash>
            scalarTable;


private:

    // Private data

        //- Phase 1
        const myPhaseModel& phase1_;

        //- Phase 2
        const myPhaseModel& phase2_;

        //- Gravitational acceleration
        const dimensionedVector& g_;

        //- Surface tension coefficient
        const dimensionedScalar sigma_;


    // Private member functions

        // Etvos number for given diameter
        tmp<volScalarField> EoH(const volScalarField& d) const;

public:

    // Constructors

        //- Construct from two phases, gravity and surface tension table
        myPhasePair
        (
            const myPhaseModel& phase1,
            const myPhaseModel& phase2,
            const dimensionedVector& g,
            const scalarTable& sigmaTable,
            const bool ordered = false
        );


    //- Destructor
    virtual ~myPhasePair();


    // Member Functions

        //- Dispersed phase
        virtual const myPhaseModel& dispersed() const;

        //- Continuous phase
        virtual const myPhaseModel& continuous() const;

        //- Pair name
        virtual word name() const;

        //- Average density
        tmp<volScalarField> rho() const;

        //- Relative velocity magnitude
        tmp<volScalarField> magUr() const;

        //- Relative velocity
        tmp<volVectorField> Ur() const;

        //- Reynolds number
        tmp<volScalarField> Re() const;

        //- Prandtl number
        tmp<volScalarField> Pr() const;

        //- Eotvos number
        tmp<volScalarField> Eo() const;

        //- Eotvos number based on hydraulic diameter type 1
        tmp<volScalarField> EoH1() const;

        //- Eotvos number based on hydraulic diameter type 2
        tmp<volScalarField> EoH2() const;

        //- Morton Number
        tmp<volScalarField> Mo() const;

        //- Takahashi Number
        tmp<volScalarField> Ta() const;

        //- Aspect ratio
        virtual tmp<volScalarField> E() const;

        // Access

            // Phase 1
            inline const myPhaseModel& phase1() const;

            // Phase 2
            inline const myPhaseModel& phase2() const;

            // Gravitational acceleration
            inline const dimensionedVector& g() const;

            // Surface tension coefficient
            inline const dimensionedScalar& sigma() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "myPhasePairI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
