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
    Foam::myPhasePairKey

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef myPhasePairKey_H
#define myPhasePairKey_H

#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class myPhasePairKey;

bool operator==(const myPhasePairKey& a, const myPhasePairKey& b);
bool operator!=(const myPhasePairKey& a, const myPhasePairKey& b);

Istream& operator>>(Istream& is, myPhasePairKey& key);
Ostream& operator<<(Ostream& os, const myPhasePairKey& key);

/*---------------------------------------------------------------------------*\
                         Class myPhasePairKey Declaration
\*---------------------------------------------------------------------------*/

class myPhasePairKey
:
    public Pair<word>
{
public:

        class hash
        :
            public Hash<myPhasePairKey>
        {
        public:

            // Constructors

                // Construct null
                hash();


            // Member operators

                // Generate a hash from a phase pair key
                label operator()(const myPhasePairKey& key) const;
        };


private:

    // Private data

        //- Flag to indicate whether ordering is important
        bool ordered_;

public:


    // Constructors

        //- Construct null
        myPhasePairKey();

        //- Construct from names and the ordering flag
        myPhasePairKey
        (
            const word& name1,
            const word& name2,
            const bool ordered
        );


    //- Destructor
    virtual ~myPhasePairKey();


    // Access

        //- Return the ordered flag
        bool ordered() const;


    // Friend Operators

        //- Test for equality
        friend bool operator==(const myPhasePairKey& a, const myPhasePairKey& b);

        //- Test for inequality
        friend bool operator!=(const myPhasePairKey& a, const myPhasePairKey& b);

        //- Read from Istream
        friend Istream& operator>>(Istream& is, myPhasePairKey& key);

        //- Write to Ostream
        friend Ostream& operator<<(Ostream& os, const myPhasePairKey& key);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
