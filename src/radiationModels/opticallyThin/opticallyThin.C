/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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

#include "opticallyThin.H"
#include "basicThermo.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(opticallyThin, 0);
        addToRadiationRunTimeSelectionTables(opticallyThin);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::radiation::opticallyThin::aPlanck
(
    const word& specie,
    const scalar Tcell
)
{
    // The range the fits were made for
    const scalar T = min(max(Tcell, scalar(300)), scalar(2500));
    const scalar x = 1000/T;

    if (specie == "H2O")
    {
        return
            -0.23093 + x*(-1.12390 + x*(9.41530 + x*(-2.99880
          + x*(0.51382 - x*1.86840e-5))));
    }
    else if (specie == "CO2")
    {
        return
            18.741 + x*(-121.310 + x*(273.500 + x*(-194.050
          + x*(56.310 - x*5.8169))));
    }
    else if (specie == "CO")
    {
        if (T <= 750)
        {
            return
                4.7869 + T*(-0.06953 + T*(2.95775e-4
              + T*(-4.25732e-7 + T*2.02894e-10)));
        }

        return
            10.09 + T*(-0.01183 + T*(4.7753e-6
          + T*(-5.87209e-10 - T*2.5334e-14)));
    }
    else if (specie == "CH4")
    {
        return
            6.6334 + T*(-0.0035686 + T*(1.6682e-8
          + T*(2.5611e-10 - T*2.6558e-14)));
    }

    FatalErrorInFunction
        << "No Planck mean absorption coefficient for " << specie << nl
        << "Available: H2O CO2 CO CH4" << exit(FatalError);

    return 0;
}


Foam::scalar Foam::radiation::opticallyThin::W(const word& specie)
{
    if (specie == "H2O") return 18.01528;
    if (specie == "CO2") return 44.0095;
    if (specie == "CO") return 28.0101;
    if (specie == "CH4") return 16.04246;

    FatalErrorInFunction
        << "No molar mass for " << specie << exit(FatalError);

    return 0;
}


void Foam::radiation::opticallyThin::readCoeffs()
{
    Tambient_.value() = coeffs_.getOrDefault<scalar>("Tambient", 300);

    species_ = coeffs_.getOrDefault<wordList>
    (
        "species",
        wordList({"H2O", "CO2", "CO", "CH4"})
    );

    for (const word& specie : species_)
    {
        // Stops here for a specie without a fit
        aPlanck(specie, 1000);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::opticallyThin::opticallyThin(const volScalarField& T)
:
    radiationModel(typeName, T),
    Tambient_("Tambient", dimTemperature, 300),
    species_(),
    a_
    (
        IOobject
        (
            "aPlanck",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    )
{
    readCoeffs();
}


Foam::radiation::opticallyThin::opticallyThin
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    Tambient_("Tambient", dimTemperature, 300),
    species_(),
    a_
    (
        IOobject
        (
            "aPlanck",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, Zero)
    )
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::opticallyThin::read()
{
    if (radiationModel::read())
    {
        readCoeffs();
        return true;
    }

    return false;
}


void Foam::radiation::opticallyThin::calculate()
{
    const auto& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    const scalarField& p = thermo.p().primitiveField();
    const scalarField& T = T_.primitiveField();

    // Mean molar mass [kg/kmol]
    const volScalarField Wmix(thermo.W());

    scalarField& a = a_.primitiveFieldRef();
    a = 0;

    for (const word& specie : species_)
    {
        const auto* Yptr = mesh_.cfindObject<volScalarField>(specie);

        if (!Yptr)
        {
            continue;
        }

        const scalarField& Y = Yptr->primitiveField();
        const scalar Wi = W(specie);

        forAll(a, celli)
        {
            // Partial pressure [atm]
            const scalar pi =
                Y[celli]*Wmix[celli]/Wi*p[celli]/101325.0;

            a[celli] += pi*aPlanck(specie, T[celli]);
        }
    }

    a_.correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField> Foam::radiation::opticallyThin::Rp() const
{
    return volScalarField::New
    (
        "Rp",
        IOobject::NO_REGISTER,
        4.0*physicoChemical::sigma*a_
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::radiation::opticallyThin::Ru() const
{
    return volScalarField::Internal::New
    (
        "Ru",
        IOobject::NO_REGISTER,
        4.0*physicoChemical::sigma*a_()*pow4(Tambient_)
    );
}


// ************************************************************************* //
