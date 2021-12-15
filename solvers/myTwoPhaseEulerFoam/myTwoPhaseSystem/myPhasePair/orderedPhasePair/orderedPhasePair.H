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
    Foam::orderedPhasePair

Description

SourceFiles
    orderedPhasePair.C

\*---------------------------------------------------------------------------*/

#ifndef orderedPhasePair_H
#define orderedPhasePair_H

#include "myPhasePair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myAspectRatioModel;

/*---------------------------------------------------------------------------*\
                         Class orderedPhasePair Declaration
\*---------------------------------------------------------------------------*/

class orderedPhasePair
:
    public myPhasePair
{
    // Private data

        //- Aspect ratio model
        autoPtr<myAspectRatioModel> aspectRatio_;


public:

    // Constructors

        //- Construct from two phases, gravity, surface tension and aspect
        //  ratio tables
        orderedPhasePair
        (
            const myPhaseModel& dispersed,
            const myPhaseModel& continuous,
            const dimensionedVector& g,
            const scalarTable& sigmaTable,
            const dictTable& aspectRatioTable
        );


    //- Destructor
    virtual ~orderedPhasePair();


    // Member Functions

        //- Dispersed phase
        virtual const myPhaseModel& dispersed() const;

        //- Continuous phase
        virtual const myPhaseModel& continuous() const;

        //- Pair name
        virtual word name() const;

        //- Aspect ratio
        virtual tmp<volScalarField> E() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
