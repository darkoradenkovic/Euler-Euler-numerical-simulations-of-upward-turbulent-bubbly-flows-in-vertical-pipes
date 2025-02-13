/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2014-2015 OpenFOAM Foundation
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
    Foam::myWallLubricationModels::Antal

Description
    Wall lubrication model of Antal et al.

    Reference:
    \verbatim
        "Analysis of phase distribution in fully developed laminar bubbly
         two-phase flow"
        Antal, S.P., Lahey Jr, R.T., and Flaherty, J.E.
        International Journal of Multiphase Flow
        Volume 17, Issue 5, September 1991, pp. 635-652
    \endverbatim

SourceFiles
    Antal.C

\*---------------------------------------------------------------------------*/

#ifndef Antal_H
#define Antal_H

#include "myWallLubricationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myPhasePair;

namespace myWallLubricationModels
{

/*---------------------------------------------------------------------------*\
                            Class Antal Declaration
\*---------------------------------------------------------------------------*/

class Antal
:
    public myWallLubricationModel
{
    // Private data

        //- Coefficient 1
        const dimensionedScalar Cw1_;

        //- Coefficient 2
        const dimensionedScalar Cw2_;


public:

    //- Runtime type information
    TypeName("Antal");


    // Constructors

        //- Construct from components
        Antal
        (
            const dictionary& dict,
            const myPhasePair& pair
        );


    //- Destructor
    virtual ~Antal();


    // Member Functions

        //- Return phase-intensive wall lubrication force
        tmp<volVectorField> Fi() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace myWallLubricationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
