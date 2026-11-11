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

Class
    Foam::radiation::opticallyThin

Group
    grpRadiationModels

Description
    Optically thin radiation: the gas emits and nothing it emits is absorbed
    again,

        Sh = - 4 sigma a (T^4 - Tambient^4)
        a  = p (X_H2O a_H2O + X_CO2 a_CO2 + X_CO a_CO + X_CH4 a_CH4)

    with the Planck mean absorption coefficients a_i(T) [1/m/atm] of the TNF
    workshop radiation model (Barlow et al., Combust. Flame 127 (2001) 2102;
    curve fits to RADCAL, valid between 300 and 2500 K). Species the
    mechanism does not have are left out.

Usage
    In constant/radiationProperties:

    \verbatim
    radiationModel  opticallyThin;

    opticallyThinCoeffs
    {
        Tambient    293;                    // default 300
        species     (H2O CO2 CO CH4);       // default these four
    }
    \endverbatim

    and in system/controlDict

    \verbatim
    libs (libopticallyThinRadiation.so);
    \endverbatim

SourceFiles
    opticallyThin.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_radiation_opticallyThin_H
#define Foam_radiation_opticallyThin_H

#include "radiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                        Class opticallyThin Declaration
\*---------------------------------------------------------------------------*/

class opticallyThin
:
    public radiationModel
{
    // Private Data

        //- Temperature of the surroundings [K]
        dimensionedScalar Tambient_;

        //- Radiating species
        wordList species_;

        //- Planck mean absorption coefficient of the mixture [1/m]
        volScalarField a_;


    // Private Member Functions

        //- Planck mean absorption coefficient of a specie [1/m/atm]
        static scalar aPlanck(const word& specie, const scalar T);

        //- Molar mass of a radiating specie [kg/kmol]
        static scalar W(const word& specie);

        //- Read the coefficients
        void readCoeffs();

        //- No copy construct
        opticallyThin(const opticallyThin&) = delete;

        //- No copy assignment
        void operator=(const opticallyThin&) = delete;


public:

    //- Runtime type information
    TypeName("opticallyThin");


    // Constructors

        //- Construct from components
        opticallyThin(const volScalarField& T);

        //- Construct from components
        opticallyThin(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~opticallyThin() = default;


    // Member Functions

        //- Update the absorption coefficient
        void calculate();

        //- Read radiationProperties dictionary
        bool read();

        //- Source term component (for power of T^4)
        tmp<volScalarField> Rp() const;

        //- Source term component (constant)
        tmp<volScalarField::Internal> Ru() const;

        //- Number of bands
        label nBands() const
        {
            return 1;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
