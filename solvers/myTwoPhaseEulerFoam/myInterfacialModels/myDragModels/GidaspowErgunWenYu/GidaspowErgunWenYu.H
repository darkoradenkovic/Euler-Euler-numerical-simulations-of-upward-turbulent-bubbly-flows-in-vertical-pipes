/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::myDragModels::GidaspowErgunWenYu

Description
    Gidaspow, Ergun, Wen and Yu drag model

    Reference:
    \verbatim
        "Multiphase flow and fluidization",
        Gidaspow, D.,
        Academic Press, New York, 1994.
    \endverbatim

SourceFiles
    GidaspowErgunWenYu.C

\*---------------------------------------------------------------------------*/

#ifndef GidaspowErgunWenYu_H
#define GidaspowErgunWenYu_H

#include "myDragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myPhasePair;

namespace myDragModels
{

class Ergun;
class WenYu;

/*---------------------------------------------------------------------------*\
                           Class GidaspowErgunWenYu Declaration
\*---------------------------------------------------------------------------*/

class GidaspowErgunWenYu
:
    public myDragModel
{
    // Private data

        //- Ergun drag model
        autoPtr<Ergun> Ergun_;

        //- Wen Yu drag model
        autoPtr<WenYu> WenYu_;


public:

    //- Runtime type information
    TypeName("GidaspowErgunWenYu");


    // Constructors

        //- Construct from a dictionary and an ordered phase pair
        GidaspowErgunWenYu
        (
            const dictionary& interfaceDict,
            const myPhasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~GidaspowErgunWenYu();


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
