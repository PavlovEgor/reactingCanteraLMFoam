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

#include "massFluxInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massFluxInletVelocityFvPatchVectorField::
massFluxInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    massFlux_(nullptr),
    rhoName_("rho"),
    rhoInlet_(-VGREAT)
{}


Foam::massFluxInletVelocityFvPatchVectorField::
massFluxInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, IOobjectOption::NO_READ),
    massFlux_(Function1<scalar>::New("massFlux", dict, &db())),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    rhoInlet_(dict.getOrDefault<scalar>("rhoInlet", -VGREAT))
{
    if (!this->readValueEntry(dict))
    {
        fvPatchVectorField::operator=(vector::zero);
    }
}


Foam::massFluxInletVelocityFvPatchVectorField::
massFluxInletVelocityFvPatchVectorField
(
    const massFluxInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    massFlux_(ptf.massFlux_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_)
{}


Foam::massFluxInletVelocityFvPatchVectorField::
massFluxInletVelocityFvPatchVectorField
(
    const massFluxInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    massFlux_(ptf.massFlux_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_)
{}


Foam::massFluxInletVelocityFvPatchVectorField::
massFluxInletVelocityFvPatchVectorField
(
    const massFluxInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    massFlux_(ptf.massFlux_.clone()),
    rhoName_(ptf.rhoName_),
    rhoInlet_(ptf.rhoInlet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::massFluxInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar m0 = massFlux_->value(db().time().timeOutputValue());

    const vectorField n(patch().nf());

    if (db().foundObject<volScalarField>(rhoName_))
    {
        const auto& rhop = patch().lookupPatchField<volScalarField>(rhoName_);

        operator==(-n*m0/rhop);
    }
    else if (rhoInlet_ > 0)
    {
        operator==(-n*m0/rhoInlet_);
    }

    // Otherwise the density does not exist yet: keep the value read

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::massFluxInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    massFlux_->writeData(os);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<scalar>("rhoInlet", -VGREAT, rhoInlet_);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        massFluxInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
