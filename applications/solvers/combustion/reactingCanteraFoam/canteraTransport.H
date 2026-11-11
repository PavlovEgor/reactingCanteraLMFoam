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
    Foam::canteraTransport

Description
    Detailed molecular transport for reactingCanteraFoam, with the transport
    coefficients from Cantera.

    The diffusive mass flux of specie k is

        j_k = - rho D_k (grad(Y_k) + Y_k grad(W)/W) - DT_k grad(T)/T + Y_k jc

    the Hirschfelder-Curtiss (mixture-averaged) approximation, written for
    the gradient of the mole fraction, plus thermal diffusion (Soret) and the
    correction flux jc = - sum_k j_k that makes the fluxes sum to zero. With
    \c diffusion \c multicomponent the first term is instead

        j_k = rho W_k/W^2 sum_j W_j D_kj grad(X_j)

    with the multicomponent diffusion coefficients D_kj, evaluated on the
    faces.

    For the specie equations the flux is split into Fick's part, which is
    implicit, and the rest, which is not:

        div(j_k) = - laplacian(rhoD_k, Y_k) + div(phiD_k) + div(phiDc Y_k)

    and the energy equation gets the heat carried by the fluxes, hsFlux().
    The fields rhoD_<specie>, phiD_<specie> and phiDc are registered, for
    boundary conditions that need the diffusive flux (burnerSpecieFlux).

    Thermodynamics and chemistry stay OpenFOAM's; the Cantera mechanism has
    to have the species of the case, by name.

Usage
    constant/canteraProperties:

    \verbatim
    mechanism   "<constant>/grimech/gri30.yaml";
    phase       gri30;                  // default: the first in the file
    diffusion   mixtureAveraged;        // unityLewis | mixtureAveraged
                                        // | multicomponent
    Soret       on;
    SoretInterval 5;                    // default 1
    \endverbatim

    \table
        diffusion       | D_k                    | Cost per step
        unityLewis      | kappa/(rho Cp)         | one evaluation per cell
        mixtureAveraged | Hirschfelder-Curtiss   | the same
        multicomponent  | D_kj, on the faces     | an N x N inverse per face
    \endtable

    Soret needs the multicomponent thermal diffusion coefficients, a
    3N x 3N system per cell (per face with multicomponent diffusion). With
    mixture-averaged diffusion they are kept for SoretInterval time steps.

SourceFiles
    canteraTransport.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_canteraTransport_H
#define Foam_canteraTransport_H

#include "rhoReactionThermo.H"
#include "fvMatrices.H"
#include "Enum.H"
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class canteraTransport Declaration
\*---------------------------------------------------------------------------*/

class canteraTransport
{
public:

    //- Diffusion models
    enum class diffusionType
    {
        unityLewis,
        mixtureAveraged,
        multicomponent
    };

    //- Names of the diffusion models
    static const Enum<diffusionType> diffusionTypeNames;


private:

    // Private Data

        //- The Cantera objects; kept out of this header, as Cantera's and
        //- OpenFOAM's macros do not go together
        struct canteraData;

        const fvMesh& mesh_;

        const rhoReactionThermo& thermo_;

        //- Mass fractions
        const PtrList<volScalarField>& Y_;

        //- constant/canteraProperties
        IOdictionary dict_;

        diffusionType diffusion_;

        bool Soret_;

        //- Time steps between updates of the thermal diffusion coefficients
        label SoretInterval_;

        //- Time index of their last update
        label SoretTimeIndex_;

        std::unique_ptr<canteraData> ct_;

        //- Dynamic viscosity [kg/m/s]
        volScalarField mu_;

        //- Thermal conductivity [W/m/K]
        volScalarField kappa_;

        //- kappa/Cp [kg/m/s]
        volScalarField alpha_;

        //- rho D_k [kg/m/s]
        PtrList<volScalarField> rhoD_;

        //- Thermal diffusion coefficients [kg/m/s]; mixture-averaged
        //- diffusion with Soret only
        PtrList<volScalarField> DT_;

        //- Explicit part of the diffusive flux of each specie [kg/s]
        PtrList<surfaceScalarField> phiD_;

        //- Correction flux [kg/s]
        surfaceScalarField phiDc_;


    // Private Member Functions

        //- Multicomponent fluxes, on the faces
        void multicomponentFluxes(PtrList<surfaceScalarField>& J) const;


public:

    // Constructors

        //- Construct from the thermo
        explicit canteraTransport(rhoReactionThermo& thermo);

        //- No copy construct
        canteraTransport(const canteraTransport&) = delete;

        //- No copy assignment
        void operator=(const canteraTransport&) = delete;


    //- Destructor
    ~canteraTransport();


    // Member Functions

        //- Dynamic viscosity [kg/m/s]
        const volScalarField& mu() const noexcept
        {
            return mu_;
        }

        //- Thermal conductivity [W/m/K]
        const volScalarField& kappa() const noexcept
        {
            return kappa_;
        }

        //- kappa/Cp [kg/m/s]
        const volScalarField& alpha() const noexcept
        {
            return alpha_;
        }

        //- rho D of a specie [kg/m/s]
        const volScalarField& rhoD(const label speciei) const
        {
            return rhoD_[speciei];
        }

        //- Transport coefficients for the current state
        void correct();

        //- Explicit and correction fluxes for the current fields
        void correctFluxes();

        //- Divergence of the diffusive flux of a specie
        tmp<fvScalarMatrix> divj(const label speciei, volScalarField& Yi) const;

        //- Sensible enthalpy carried by the diffusive fluxes, less what the
        //- enthalpy Laplacian with kappa/Cp already carries [W]
        tmp<surfaceScalarField> hsFlux() const;

        //- Viscous term of the momentum equation
        tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
