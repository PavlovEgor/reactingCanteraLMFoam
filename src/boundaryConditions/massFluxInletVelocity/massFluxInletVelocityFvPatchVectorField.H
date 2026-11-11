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
    Foam::massFluxInletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Inlet velocity from a mass flux per unit area,

        U_b = - n m0/rho_b

    so the mass entering through the patch stays m0 whatever the density on
    the patch becomes, e.g. when a flame heats the surface of a burner.

Usage
    \table
        Property    | Description                         | Required | Default
        massFlux    | mass flux [kg/m2/s], Function1 of time | yes   |
        rho         | density field name                  | no       | rho
        rhoInlet    | density used until rho exists       | no       | none
    \endtable

    \verbatim
    <patchName>
    {
        type        massFluxInletVelocity;
        massFlux    0.1;
        value       uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    massFluxInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef massFluxInletVelocityFvPatchVectorField_H
#define massFluxInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class massFluxInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class massFluxInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Mass flux per unit area [kg/m2/s]
        autoPtr<Function1<scalar>> massFlux_;

        //- Name of the density field
        word rhoName_;

        //- Density to use while the density field does not exist yet
        scalar rhoInlet_;


public:

    //- Runtime type information
    TypeName("massFluxInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        massFluxInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        massFluxInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        massFluxInletVelocityFvPatchVectorField
        (
            const massFluxInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        massFluxInletVelocityFvPatchVectorField
        (
            const massFluxInletVelocityFvPatchVectorField&
        );

        //- Construct as copy setting internal field reference
        massFluxInletVelocityFvPatchVectorField
        (
            const massFluxInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new massFluxInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new massFluxInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

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
