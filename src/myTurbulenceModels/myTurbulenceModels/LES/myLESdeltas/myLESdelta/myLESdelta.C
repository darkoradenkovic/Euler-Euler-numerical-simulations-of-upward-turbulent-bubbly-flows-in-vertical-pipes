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

\*---------------------------------------------------------------------------*/

#include "myLESdelta.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(myLESdelta, 0);
    defineRunTimeSelectionTable(myLESdelta, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::myLESdelta::myLESdelta
(
    const word& name,
    const myTurbulenceModel& turbulence
)
:
    myTurbulenceModel_(turbulence),
    delta_
    (
        IOobject
        (
            name,
            turbulence.mesh().time().timeName(),
            turbulence.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        turbulence.mesh(),
        dimensionedScalar(name, dimLength, SMALL),
        calculatedFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::myLESdelta> Foam::myLESdelta::New
(
    const word& name,
    const myTurbulenceModel& turbulence,
    const dictionary& dict,
    const word& lookupName
)
{
    const word deltaType(dict.get<word>(lookupName));

    Info<< "Selecting LES " << lookupName << " type " << deltaType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(deltaType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown myLESdelta type "
            << deltaType << nl << nl
            << "Valid myLESdelta types :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<myLESdelta>(cstrIter()(name, turbulence, dict));
}


Foam::autoPtr<Foam::myLESdelta> Foam::myLESdelta::New
(
    const word& name,
    const myTurbulenceModel& turbulence,
    const dictionary& dict,
    const dictionaryConstructorTable& additionalConstructors,
    const word& lookupName
)
{
    const word deltaType(dict.get<word>(lookupName));

    Info<< "Selecting LES " << lookupName << " type " << deltaType << endl;

    // First any additional ones
    {
        auto cstrIter = additionalConstructors.cfind(deltaType);

        if (cstrIter.found())
        {
            return autoPtr<myLESdelta>(cstrIter()(name, turbulence, dict));
        }
    }

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(deltaType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown myLESdelta type "
            << deltaType << nl << nl
            << "Valid myLESdelta types :" << endl
            << additionalConstructors.sortedToc()
            << " and "
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<myLESdelta>(cstrIter()(name, turbulence, dict));
}


// ************************************************************************* //
