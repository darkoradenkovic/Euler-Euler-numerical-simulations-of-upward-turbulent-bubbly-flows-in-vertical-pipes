/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::myLESdelta

Description
    Abstract base class for LES deltas

SourceFiles
    myLESdelta.C

\*---------------------------------------------------------------------------*/

#ifndef myLESdelta_H
#define myLESdelta_H

#include "myTurbulenceModel.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class myLESdelta Declaration
\*---------------------------------------------------------------------------*/

class myLESdelta
{

protected:

    // Protected data

        const myTurbulenceModel& myTurbulenceModel_;

        volScalarField delta_;


    // Protected Member Functions

        //- No copy construct
        myLESdelta(const myLESdelta&) = delete;

        //- No copy assignment
        void operator=(const myLESdelta&) = delete;


public:

    //- Runtime type information
    TypeName("myLESdelta");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            myLESdelta,
            dictionary,
            (
                const word& name,
                const myTurbulenceModel& turbulence,
                const dictionary& dict
            ),
            (name, turbulence, dict)
        );


    // Constructors

        //- Construct from name, myTurbulenceModel and dictionary
        myLESdelta
        (
            const word& name,
            const myTurbulenceModel& turbulence
        );


    // Selectors

        //- Return a reference to the selected LES delta
        static autoPtr<myLESdelta> New
        (
            const word& name,
            const myTurbulenceModel& turbulence,
            const dictionary& dict,
            const word& lookupName = "delta"
        );

        //- Return a reference to the selected LES delta
        static autoPtr<myLESdelta> New
        (
            const word& name,
            const myTurbulenceModel& turbulence,
            const dictionary& dict,
            const dictionaryConstructorTable& additionalConstructors,
            const word& lookupName = "delta"
        );


    //- Destructor
    virtual ~myLESdelta()
    {}


    // Member Functions

        //- Return myTurbulenceModel reference
        const myTurbulenceModel& turbulence() const
        {
            return myTurbulenceModel_;
        }

        //- Read the myLESdelta dictionary
        virtual void read(const dictionary&) = 0;

        // Correct values
        virtual void correct() = 0;


    // Member Operators

        virtual operator const volScalarField&() const
        {
            return delta_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
