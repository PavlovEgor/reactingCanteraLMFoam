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

#include "burnerSpecieFluxFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::burnerSpecieFluxFvPatchScalarField::burnerSpecieFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    Y0_(0),
    phiName_("phi")
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::burnerSpecieFluxFvPatchScalarField::burnerSpecieFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    Y0_(dict.get<scalar>("Y0")),
    phiName_(dict.getOrDefault<word>("phi", "phi"))
{
    fvPatchFieldBase::readDict(dict);

    refValue() = Y0_;
    refGrad() = 0;
    valueFraction() = 1;

    if (!this->readValueEntry(dict))
    {
        fvPatchScalarField::operator=(Y0_);
    }
}


Foam::burnerSpecieFluxFvPatchScalarField::burnerSpecieFluxFvPatchScalarField
(
    const burnerSpecieFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    Y0_(ptf.Y0_),
    phiName_(ptf.phiName_)
{}


Foam::burnerSpecieFluxFvPatchScalarField::burnerSpecieFluxFvPatchScalarField
(
    const burnerSpecieFluxFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    Y0_(ptf.Y0_),
    phiName_(ptf.phiName_)
{}


Foam::burnerSpecieFluxFvPatchScalarField::burnerSpecieFluxFvPatchScalarField
(
    const burnerSpecieFluxFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    Y0_(ptf.Y0_),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::burnerSpecieFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const word& specieName = internalField().name();
    const word rhoDName("rhoD_" + specieName);
    const word phiDName("phiD_" + specieName);

    refGrad() = 0;

    if
    (
        !db().foundObject<surfaceScalarField>(phiName_)
     || !db().foundObject<volScalarField>(rhoDName)
     || !db().foundObject<surfaceScalarField>(phiDName)
     || !db().foundObject<surfaceScalarField>("phiDc")
    )
    {
        // No flux or no diffusion model yet: the supplied mixture
        refValue() = Y0_;
        valueFraction() = 1;

        mixedFvPatchScalarField::updateCoeffs();
        return;
    }

    const auto& phip = patch().lookupPatchField<surfaceScalarField>(phiName_);
    const auto& rhoDp = patch().lookupPatchField<volScalarField>(rhoDName);
    const auto& phiDp = patch().lookupPatchField<surfaceScalarField>(phiDName);
    const auto& phiDcp = patch().lookupPatchField<surfaceScalarField>("phiDc");

    const scalarField& magSf = patch().magSf();
    const scalarField& delta = patch().deltaCoeffs();

    scalarField& ref = refValue();
    scalarField& f = valueFraction();

    forAll(f, facei)
    {
        // Mass flux into the domain
        const scalar m0 = -phip[facei]/magSf[facei];

        if (m0 <= SMALL)
        {
            // Outflow
            ref[facei] = Y0_;
            f[facei] = 0;
            continue;
        }

        // The correction flux carries the specie as the flow does; it is a
        // small part of m0 unless the flame is about to enter the burner
        const scalar m = max(m0 - phiDcp[facei]/magSf[facei], 0.1*m0);
        const scalar e = phiDp[facei]/magSf[facei];

        ref[facei] = min(max((m0*Y0_ + e)/m, scalar(0)), scalar(1));
        f[facei] = m/(m + rhoDp[facei]*delta[facei]);
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::burnerSpecieFluxFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("Y0", Y0_);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        burnerSpecieFluxFvPatchScalarField
    );
}

// ************************************************************************* //
