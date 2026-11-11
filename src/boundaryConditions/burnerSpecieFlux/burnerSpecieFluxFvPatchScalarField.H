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
    Foam::burnerSpecieFluxFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    Robin condition for a specie mass fraction on the surface of a burner:
    the total flux of the specie, convection plus diffusion, is that of the
    unburnt mixture,

        m0 (Y_b - Y0) n = - j.n

    with m0 the mass flux through the patch, Y0 the mass fraction in the
    supplied mixture and j the diffusive flux of the specie. A flame sitting
    close to the surface then changes the composition on the surface instead
    of the amount of each specie that enters (Margolis, Combust. Sci.
    Technol. 22 (1980); Li et al., Combust. Flame 276 (2025) 114115, Eq. 2).

    The diffusive flux is the one reactingCanteraFoam solves with,

        j.n Sf = - rhoD snGrad(Y) Sf + phiD + phiDc Y_b

    Fick's part, the explicit part (molar mass gradient, Soret,
    multicomponent) and the correction flux, from the fields \c rhoD_<specie>,
    \c phiD_<specie> and \c phiDc that the solver registers. That gives a
    mixed condition with

        valueFraction = (m0 - c)/(m0 - c + rhoD deltaCoeffs)
        refValue      = (m0 Y0 + e)/(m0 - c)

    where c = phiDc/Sf and e = phiD/Sf. Until those fields exist, and where
    the flow leaves the domain, it is a fixed value Y0 and a zero gradient,
    respectively.

Usage
    \table
        Property | Description                            | Required | Default
        Y0       | mass fraction in the supplied mixture  | yes      |
        phi      | mass flux field name                   | no       | phi
    \endtable

    \verbatim
    <patchName>
    {
        type        burnerSpecieFlux;
        Y0          0.0629;
        value       uniform 0.0629;
    }
    \endverbatim

SourceFiles
    burnerSpecieFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef burnerSpecieFluxFvPatchScalarField_H
#define burnerSpecieFluxFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class burnerSpecieFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class burnerSpecieFluxFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Mass fraction in the supplied mixture
        scalar Y0_;

        //- Name of the mass flux field
        word phiName_;


public:

    //- Runtime type information
    TypeName("burnerSpecieFlux");


    // Constructors

        //- Construct from patch and internal field
        burnerSpecieFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        burnerSpecieFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        burnerSpecieFluxFvPatchScalarField
        (
            const burnerSpecieFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        burnerSpecieFluxFvPatchScalarField
        (
            const burnerSpecieFluxFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        burnerSpecieFluxFvPatchScalarField
        (
            const burnerSpecieFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new burnerSpecieFluxFvPatchScalarField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new burnerSpecieFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Mass fraction in the supplied mixture
        scalar Y0() const noexcept
        {
            return Y0_;
        }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
