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
    Foam::myDragModels::TomiyamaCorrelatedLM

Description
    Correlation of Tomiyama et al.

    Reference:
    \verbatim
        "Terminal velocity of single bubbles in surface tension force dominant
         regime"
        Tomiyama, T., Celata, G.P., Hosokawa, S., Yoshida, S.,
        International Journal of Multiphase Flow
        Volume 28, Issue 9, September 2002, pp. 1497-1519
    \endverbatim

SourceFiles
    TomiyamaCorrelatedLM.C

\*---------------------------------------------------------------------------*/

#ifndef TomiyamaCorrelatedLM_H
#define TomiyamaCorrelatedLM_H

#include "myDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myPhasePair;

namespace myDragModels
{

/*---------------------------------------------------------------------------*\
                           Class TomiyamaCorrelatedLM Declaration
\*---------------------------------------------------------------------------*/

class TomiyamaCorrelatedLM
:
    public myDragModel
{
    // Private data

        //- Coefficient
        const dimensionedScalar A_;


public:

    //- Runtime type information
    TypeName("TomiyamaCorrelatedLM");


    // Constructors

        //- Construct from a dictionary and a phase pair
        TomiyamaCorrelatedLM
        (
            const dictionary& dict,
            const myPhasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~TomiyamaCorrelatedLM();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace myDragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
