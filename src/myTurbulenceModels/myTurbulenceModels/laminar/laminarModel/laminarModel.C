/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "laminarModel.H"
#include "Stokes.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
void Foam::laminarModel<BasicMyTurbulenceModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
Foam::laminarModel<BasicMyTurbulenceModel>::laminarModel
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicMyTurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    laminarDict_(this->subOrEmptyDict("laminar")),
    printCoeffs_(laminarDict_.lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(laminarDict_.optionalSubDict(type + "Coeffs"))
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
Foam::autoPtr<Foam::laminarModel<BasicMyTurbulenceModel>>
Foam::laminarModel<BasicMyTurbulenceModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    IOdictionary modelDict
    (
        IOobject
        (
            IOobject::groupName(propertiesName, alphaRhoPhi.group()),
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    if (modelDict.found("laminar"))
    {
        // get model name, but do not register the dictionary
        // otherwise it is registered in the database twice
        const word modelType
        (
            modelDict.subDict("laminar").get<word>("laminarModel")
        );

        Info<< "Selecting laminar stress model " << modelType << endl;

        auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

        if (!cstrIter.found())
        {
            FatalErrorInFunction
                << "Unknown laminarModel type "
                << modelType << nl << nl
                << "Valid laminarModel types :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<laminarModel>
        (
            cstrIter()
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                transport, propertiesName)
        );
    }
    else
    {
        Info<< "Selecting laminar stress model "
            << laminarModels::Stokes<BasicMyTurbulenceModel>::typeName << endl;

        return autoPtr<laminarModel>
        (
            new laminarModels::Stokes<BasicMyTurbulenceModel>
            (
                alpha,
                rho,
                U,
                alphaRhoPhi,
                phi,
                transport,
                propertiesName
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMyTurbulenceModel>
bool Foam::laminarModel<BasicMyTurbulenceModel>::read()
{
    if (BasicMyTurbulenceModel::read())
    {
        laminarDict_ <<= this->subDict("laminar");

        coeffDict_ <<= laminarDict_.optionalSubDict(type() + "Coeffs");

        return true;
    }

    return false;
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicMyTurbulenceModel>::nut() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("nut", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(dimViscosity, Zero)
    );
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::laminarModel<BasicMyTurbulenceModel>::nut
(
    const label patchi
) const
{
    return tmp<scalarField>
    (
        new scalarField(this->mesh_.boundary()[patchi].size(), Zero)
    );
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicMyTurbulenceModel>::nuEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("nuEff", this->alphaRhoPhi_.group()), this->nu()
        )
    );
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::scalarField>
Foam::laminarModel<BasicMyTurbulenceModel>::nuEff
(
    const label patchi
) const
{
    return this->nu(patchi);
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicMyTurbulenceModel>::k() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("k", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(sqr(this->U_.dimensions()), Zero)
    );
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::volScalarField>
Foam::laminarModel<BasicMyTurbulenceModel>::epsilon() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(sqr(this->U_.dimensions())/dimTime, Zero)
    );
}


template<class BasicMyTurbulenceModel>
Foam::tmp<Foam::volSymmTensorField>
Foam::laminarModel<BasicMyTurbulenceModel>::R() const
{
    return tmp<volSymmTensorField>::New
    (
        IOobject
        (
            IOobject::groupName("R", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedSymmTensor(sqr(this->U_.dimensions()), Zero)
    );
}


template<class BasicMyTurbulenceModel>
void Foam::laminarModel<BasicMyTurbulenceModel>::correct()
{
    BasicMyTurbulenceModel::correct();
}


// ************************************************************************* //
