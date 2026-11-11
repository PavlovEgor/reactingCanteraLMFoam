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

// Cantera first: its headers do not survive OpenFOAM's macros
#include "cantera/core.h"
#include "cantera/transport/TransportFactory.h"

#include "canteraTransport.H"
#include "fvCFD.H"
#include "wallFvPatch.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum<Foam::canteraTransport::diffusionType>
Foam::canteraTransport::diffusionTypeNames
({
    { diffusionType::unityLewis, "unityLewis" },
    { diffusionType::mixtureAveraged, "mixtureAveraged" },
    { diffusionType::multicomponent, "multicomponent" },
});


// * * * * * * * * * * * * * * * Local Classes * * * * * * * * * * * * * * * //

struct Foam::canteraTransport::canteraData
{
    std::shared_ptr<Cantera::Solution> solution;
    std::shared_ptr<Cantera::ThermoPhase> gas;

    //- Mixture-averaged transport
    std::shared_ptr<Cantera::Transport> mix;

    //- Multicomponent transport, if anything needs it
    std::shared_ptr<Cantera::Transport> multi;

    //- Cantera index of each OpenFOAM specie
    std::vector<size_t> map;

    //- Molar masses in OpenFOAM order [kg/kmol]
    std::vector<double> W;

    // Work arrays, Cantera order
    std::vector<double> y, d, dt, dkj;

    size_t n() const
    {
        return y.size();
    }

    void setState(const double T, const double p)
    {
        gas->setState_TPY(T, p, y.data());
    }
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::canteraTransport::canteraTransport(rhoReactionThermo& thermo)
:
    mesh_(thermo.T().mesh()),
    thermo_(thermo),
    Y_(thermo.composition().Y()),
    dict_
    (
        IOobject
        (
            "canteraProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    diffusion_(diffusionTypeNames.get("diffusion", dict_)),
    Soret_(dict_.getOrDefault<Switch>("Soret", false)),
    SoretInterval_(max(dict_.getOrDefault<label>("SoretInterval", 1), 1)),
    SoretTimeIndex_(-1),
    ct_(new canteraData),
    mu_
    (
        IOobject("mu", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimDynamicViscosity, Zero)
    ),
    kappa_
    (
        IOobject("kappa", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimPower/dimLength/dimTemperature, Zero)
    ),
    alpha_
    (
        IOobject("alpha", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimDynamicViscosity, Zero)
    ),
    rhoD_(Y_.size()),
    DT_(),
    phiD_(Y_.size()),
    phiDc_
    (
        IOobject("phiDc", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimMass/dimTime, Zero)
    )
{
    const fileName mechanism(dict_.get<fileName>("mechanism").expand());
    const word phase(dict_.getOrDefault<word>("phase", word::null));

    const bool needMulti =
        Soret_ || diffusion_ == diffusionType::multicomponent;

    Info<< "Cantera transport" << nl
        << "    mechanism : " << mechanism << nl
        << "    diffusion : " << diffusionTypeNames[diffusion_] << nl
        << "    Soret     : " << Switch(Soret_) << nl << endl;

    canteraData& ct = *ct_;

    try
    {
        ct.solution =
            Cantera::newSolution(mechanism, phase, "mixture-averaged");
        ct.gas = ct.solution->thermo();
        ct.mix = ct.solution->transport();

        if (needMulti)
        {
            ct.multi = Cantera::newTransport(ct.gas, "multicomponent");
        }
    }
    catch (const Cantera::CanteraError& err)
    {
        FatalIOErrorInFunction(dict_)
            << "Cantera: " << err.what() << exit(FatalIOError);
    }

    const size_t n = ct.gas->nSpecies();
    const speciesTable& species = thermo.composition().species();

    ct.map.resize(species.size());
    ct.W.resize(species.size());
    ct.y.assign(n, 0);
    ct.d.assign(n, 0);
    ct.dt.assign(n, 0);

    if (diffusion_ == diffusionType::multicomponent)
    {
        ct.dkj.assign(n*n, 0);
    }

    forAll(species, i)
    {
        const size_t k = ct.gas->speciesIndex(species[i]);

        if (k == Cantera::npos)
        {
            FatalIOErrorInFunction(dict_)
                << "Specie " << species[i] << " is not in " << mechanism
                << exit(FatalIOError);
        }

        ct.map[i] = k;
        ct.W[i] = ct.gas->molecularWeight(k);

        rhoD_.set
        (
            i,
            new volScalarField
            (
                IOobject("rhoD_" + species[i], mesh_.time().timeName(), mesh_),
                mesh_,
                dimensionedScalar(dimDynamicViscosity, Zero)
            )
        );

        phiD_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject("phiD_" + species[i], mesh_.time().timeName(), mesh_),
                mesh_,
                dimensionedScalar(dimMass/dimTime, Zero)
            )
        );
    }

    if (Soret_ && diffusion_ != diffusionType::multicomponent)
    {
        DT_.resize(species.size());

        forAll(species, i)
        {
            DT_.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "DT_" + species[i],
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar(dimDynamicViscosity, Zero)
                )
            );
        }
    }

    phiDc_.setOriented();
    forAll(phiD_, i)
    {
        phiD_[i].setOriented();
    }

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::canteraTransport::~canteraTransport()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::canteraTransport::correct()
{
    canteraData& ct = *ct_;

    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    const bool multi = diffusion_ == diffusionType::multicomponent;
    const label timeIndex = mesh_.time().timeIndex();

    const bool cellDT =
        DT_.size()
     && (
            SoretTimeIndex_ < 0
         || timeIndex - SoretTimeIndex_ >= SoretInterval_
        );

    if (cellDT)
    {
        SoretTimeIndex_ = timeIndex;
    }

    // One state: the transport coefficients into element i of the fields
    auto evaluate = [&]
    (
        const scalar Ti,
        const scalar pi,
        const label i,
        scalarField& mu,
        scalarField& kappa,
        const UList<scalarField*>& rhoD,
        const UList<scalarField*>& DT
    )
    {
        ct.setState(Ti, pi);

        mu[i] = ct.mix->viscosity();
        kappa[i] =
            multi
          ? ct.multi->thermalConductivity()
          : ct.mix->thermalConductivity();

        if (diffusion_ != diffusionType::unityLewis)
        {
            ct.mix->getMixDiffCoeffs(ct.d.data());
            const scalar rho = ct.gas->density();

            forAll(rhoD, k)
            {
                (*rhoD[k])[i] = rho*ct.d[ct.map[k]];
            }
        }

        if (cellDT)
        {
            ct.multi->getThermalDiffCoeffs(ct.dt.data());

            forAll(DT, k)
            {
                (*DT[k])[i] = ct.dt[ct.map[k]];
            }
        }
    };

    List<scalarField*> rhoD(Y_.size());
    List<scalarField*> DT(DT_.size());

    try
    {
        // Cells
        {
            forAll(rhoD, k)
            {
                rhoD[k] = &rhoD_[k].primitiveFieldRef();
            }
            forAll(DT, k)
            {
                DT[k] = &DT_[k].primitiveFieldRef();
            }

            scalarField& mu = mu_.primitiveFieldRef();
            scalarField& kappa = kappa_.primitiveFieldRef();

            forAll(mu, celli)
            {
                forAll(Y_, k)
                {
                    ct.y[ct.map[k]] = max(Y_[k][celli], scalar(0));
                }

                evaluate(T[celli], p[celli], celli, mu, kappa, rhoD, DT);
            }
        }

        // Boundary faces; the coupled and constraint patches evaluate
        // themselves below
        forAll(mu_.boundaryField(), patchi)
        {
            if (!isA<calculatedFvPatchScalarField>(mu_.boundaryField()[patchi]))
            {
                continue;
            }

            forAll(rhoD, k)
            {
                rhoD[k] = &rhoD_[k].boundaryFieldRef()[patchi];
            }
            forAll(DT, k)
            {
                DT[k] = &DT_[k].boundaryFieldRef()[patchi];
            }

            scalarField& mu = mu_.boundaryFieldRef()[patchi];
            scalarField& kappa = kappa_.boundaryFieldRef()[patchi];

            const scalarField& Tp = T.boundaryField()[patchi];
            const scalarField& pp = p.boundaryField()[patchi];

            forAll(mu, facei)
            {
                forAll(Y_, k)
                {
                    ct.y[ct.map[k]] =
                        max(Y_[k].boundaryField()[patchi][facei], scalar(0));
                }

                evaluate(Tp[facei], pp[facei], facei, mu, kappa, rhoD, DT);
            }
        }
    }
    catch (const Cantera::CanteraError& err)
    {
        FatalErrorInFunction
            << "Cantera: " << err.what() << exit(FatalError);
    }

    mu_.correctBoundaryConditions();
    kappa_.correctBoundaryConditions();

    alpha_ = kappa_/thermo_.Cp();

    forAll(rhoD_, k)
    {
        if (diffusion_ == diffusionType::unityLewis)
        {
            rhoD_[k] = alpha_;
        }
        else
        {
            rhoD_[k].correctBoundaryConditions();
        }
    }

    forAll(DT_, k)
    {
        DT_[k].correctBoundaryConditions();
    }
}


void Foam::canteraTransport::multicomponentFluxes
(
    PtrList<surfaceScalarField>& J
) const
{
    canteraData& ct = *ct_;

    const volScalarField& T = thermo_.T();
    const label nSpecie = Y_.size();
    const size_t n = ct.n();

    const surfaceScalarField Tf(fvc::interpolate(T));
    const surfaceScalarField pf(fvc::interpolate(thermo_.p()));
    const surfaceScalarField snGradT(fvc::snGrad(T));

    const volScalarField W(thermo_.W());

    PtrList<surfaceScalarField> Yf(nSpecie);
    PtrList<surfaceScalarField> snGradX(nSpecie);

    forAll(Y_, k)
    {
        Yf.set(k, fvc::interpolate(Y_[k]));
        snGradX.set(k, fvc::snGrad(volScalarField(Y_[k]*W/ct.W[k])));
    }

    // One face: the fluxes per unit area into j
    scalarField j(nSpecie);

    auto evaluate = [&]
    (
        const scalar Ti,
        const scalar pi,
        const scalar snGradTi,
        const scalarField& gradX
    )
    {
        ct.setState(Ti, pi);

        ct.multi->getMultiDiffCoeffs(n, ct.dkj.data());

        if (Soret_)
        {
            ct.multi->getThermalDiffCoeffs(ct.dt.data());
        }

        const scalar Wm = ct.gas->meanMolecularWeight();
        const scalar rhoByW2 = ct.gas->density()/sqr(Wm);

        forAll(j, k)
        {
            const size_t kc = ct.map[k];

            scalar sum = 0;
            forAll(j, l)
            {
                // D_kl is element (k, l) of a column-major matrix
                sum += ct.W[l]*ct.dkj[n*ct.map[l] + kc]*gradX[l];
            }

            j[k] = rhoByW2*ct.W[k]*sum;

            if (Soret_)
            {
                j[k] -= ct.dt[kc]*snGradTi/Ti;
            }
        }
    };

    scalarField gradX(nSpecie);

    try
    {
        forAll(Tf, facei)
        {
            forAll(Y_, k)
            {
                ct.y[ct.map[k]] = max(Yf[k][facei], scalar(0));
                gradX[k] = snGradX[k][facei];
            }

            evaluate(Tf[facei], pf[facei], snGradT[facei], gradX);

            forAll(J, k)
            {
                J[k][facei] = j[k]*mesh_.magSf()[facei];
            }
        }

        forAll(Tf.boundaryField(), patchi)
        {
            const fvPatch& patch = mesh_.boundary()[patchi];

            // Coupled patches and those something can diffuse through
            if
            (
                !patch.coupled()
             && !isA<calculatedFvPatchScalarField>
                (
                    mu_.boundaryField()[patchi]
                )
            )
            {
                continue;
            }

            const scalarField& magSf = patch.magSf();

            forAll(magSf, facei)
            {
                forAll(Y_, k)
                {
                    ct.y[ct.map[k]] =
                        max(Yf[k].boundaryField()[patchi][facei], scalar(0));
                    gradX[k] = snGradX[k].boundaryField()[patchi][facei];
                }

                evaluate
                (
                    Tf.boundaryField()[patchi][facei],
                    pf.boundaryField()[patchi][facei],
                    snGradT.boundaryField()[patchi][facei],
                    gradX
                );

                forAll(J, k)
                {
                    J[k].boundaryFieldRef()[patchi][facei] =
                        j[k]*magSf[facei];
                }
            }
        }
    }
    catch (const Cantera::CanteraError& err)
    {
        FatalErrorInFunction
            << "Cantera: " << err.what() << exit(FatalError);
    }
}


void Foam::canteraTransport::correctFluxes()
{
    const volScalarField& T = thermo_.T();
    const surfaceScalarField& magSf = mesh_.magSf();
    const label nSpecie = Y_.size();

    // The whole diffusive flux of each specie, before the correction
    PtrList<surfaceScalarField> J(nSpecie);

    forAll(J, k)
    {
        J.set
        (
            k,
            new surfaceScalarField
            (
                IOobject
                (
                    "J",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh_,
                dimensionedScalar(dimMass/dimTime, Zero)
            )
        );
    }

    forAll(J, k)
    {
        J[k].setOriented();
    }

    if (diffusion_ == diffusionType::multicomponent)
    {
        multicomponentFluxes(J);
    }
    else
    {
        const volScalarField W(thermo_.W());
        const surfaceScalarField snGradWByW
        (
            fvc::snGrad(W)/fvc::interpolate(W)
        );

        tmp<surfaceScalarField> tsnGradTByT;
        if (DT_.size())
        {
            tsnGradTByT = fvc::snGrad(T)/fvc::interpolate(T);
        }

        forAll(J, k)
        {
            J[k] =
               -fvc::interpolate(rhoD_[k])*magSf
               *(
                    fvc::snGrad(Y_[k])
                  + fvc::interpolate(Y_[k])*snGradWByW
                );

            if (DT_.size())
            {
                J[k] -= fvc::interpolate(DT_[k])*tsnGradTByT()*magSf;
            }
        }
    }

    phiDc_ = dimensionedScalar(phiDc_.dimensions(), Zero);

    forAll(J, k)
    {
        // Nothing diffuses through a wall
        forAll(mesh_.boundary(), patchi)
        {
            if (isA<wallFvPatch>(mesh_.boundary()[patchi]))
            {
                J[k].boundaryFieldRef()[patchi] = Zero;
            }
        }

        phiDc_ -= J[k];

        // Less Fick's part, which the Laplacian has
        phiD_[k] = J[k] + fvc::interpolate(rhoD_[k])*fvc::snGrad(Y_[k])*magSf;
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::canteraTransport::divj
(
    const label k,
    volScalarField& Yi
) const
{
    return
    (
      - fvm::laplacian(rhoD_[k], Yi, "laplacian(rhoD,Yi)")
      + fvm::div(phiDc_, Yi, "div(phiDc,Yi)")
      + fvc::div(phiD_[k])
    );
}


Foam::tmp<Foam::surfaceScalarField> Foam::canteraTransport::hsFlux() const
{
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();
    const basicSpecieMixture& composition = thermo_.composition();
    const surfaceScalarField& magSf = mesh_.magSf();

    auto tflux = surfaceScalarField::New
    (
        "hsFlux",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimPower, Zero)
    );
    surfaceScalarField& flux = tflux.ref();

    const surfaceScalarField alphaf(fvc::interpolate(alpha_));

    volScalarField hsi
    (
        IOobject
        (
            "hsi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dimEnergy/dimMass, Zero)
    );

    forAll(Y_, k)
    {
        forAll(hsi, celli)
        {
            hsi[celli] = composition.Hs(k, p[celli], T[celli]);
        }

        forAll(hsi.boundaryField(), patchi)
        {
            scalarField& hsp = hsi.boundaryFieldRef()[patchi];
            const scalarField& pp = p.boundaryField()[patchi];
            const scalarField& Tp = T.boundaryField()[patchi];

            forAll(hsp, facei)
            {
                hsp[facei] = composition.Hs(k, pp[facei], Tp[facei]);
            }
        }

        flux +=
            fvc::interpolate(hsi)
           *(
                phiD_[k]
              + fvc::flux(phiDc_, Y_[k], "div(phiDc,Yi)")
              - (fvc::interpolate(rhoD_[k]) - alphaf)*fvc::snGrad(Y_[k])*magSf
            );
    }

    return tflux;
}


Foam::tmp<Foam::fvVectorMatrix> Foam::canteraTransport::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvc::div(mu_*dev2(T(fvc::grad(U))), "div(mu*dev2(T(grad(U))))")
      - fvm::laplacian(mu_, U, "laplacian(mu,U)")
    );
}


// ************************************************************************* //
