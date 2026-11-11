# reactingCanteraLMFoam

`reactingCanteraLMFoam`: an OpenFOAM v2412 solver for laminar flames that need detailed
molecular transport, with the boundary conditions of a flat porous burner and an optically
thin radiation model. Laid out like an OpenFOAM plugin; it has nothing to do with any
particular chemistry model, and runs with OpenFOAM's own or with
[kodes4foam](https://github.com/PavlovEgor/kodes4foam)'s `kodesChemistryModel`.

## Layout

```
Allwmake, Allwclean
src/
  boundaryConditions/           massFluxInletVelocity, burnerSpecieFlux
  radiationModels/              opticallyThin
applications/
  solvers/combustion/
    reactingCanteraLMFoam/      the solver, and canteraTransport: Cantera's transport
external/
  cantera/                      the Cantera submodule, 3.2.0
requirements.txt                SCons, for the Cantera build
```

## Building

Needs OpenFOAM v2412 or later sourced, the Cantera submodule and SCons:

```bash
git submodule update --init external/cantera
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt     # SCons
./Allwmake
```

The first `./Allwmake` builds Cantera as a static library inside `external/cantera` (a
quarter of an hour on two cores, `WM_NCOMPPROCS` sets how many) and links it into the
solver, so a case needs nothing of Cantera at run time but the mechanism file.
`CANTERA_DIR=<prefix>` uses an installed Cantera instead. SCons is `$SCONS`, the one in
`.venv`, in the active virtual environment or on the `PATH`; without it the solver is
skipped with a note saying so. The libraries go to `$FOAM_USER_LIBBIN` and the solver to
`$FOAM_USER_APPBIN`, or where `-prefix` says, as for any plugin.

## The solver

`rhoReactingBuoyantFoam` with the molecular transport a resolved laminar flame needs,
which `reactingFoam` does not have - it diffuses every specie and heat at one rate. It
follows EBIdnsFoam, the solver the flat-burner computations of
Li et al., Combust. Flame 276 (2025) 114115, were made with:

| | |
|---|---|
| Viscosity, conductivity, diffusion coefficients | Cantera, from the mechanism's YAML file |
| Diffusion | `mixtureAveraged` (Hirschfelder-Curtiss, with the correction flux), `multicomponent` (D_kj, evaluated on the faces) or `unityLewis` |
| Thermal diffusion | `Soret on`: multicomponent thermal diffusion coefficients |
| Energy | the enthalpy carried by the diffusive fluxes |
| Radiation | OpenFOAM's `radiationModel`; `opticallyThin` here is H2O, CO2, CO and CH4 with the TNF Planck mean coefficients |
| Buoyancy | `p_rgh`, as in `rhoReactingBuoyantFoam` |

Thermodynamics and chemistry stay OpenFOAM's, so it runs with `kodesChemistryModel` like
`reactingFoam` does. `constant/canteraProperties`:

```
mechanism       "<constant>/grimech/grimech.yaml";
diffusion       mixtureAveraged;    // unityLewis | mixtureAveraged | multicomponent
Soret           on;
SoretInterval   5;                  // time steps between updates of the coefficients
```

The Cantera mechanism has to have the species of the case, by name. Thermal diffusion
costs a 3N x 3N solve per cell, which is why its coefficients can be kept for a few
steps; `multicomponent` costs an N x N inverse per face on every step and is there to
check `mixtureAveraged` against, not to run with.

Two boundary conditions come with it, in `libburnerBoundaryConditions.so`:

| | |
|---|---|
| `massFluxInletVelocity` | `U = -n massFlux/rho`: an inlet set by its mass flux per area |
| `burnerSpecieFlux` | Robin condition on the surface of a burner: convection plus diffusion of a specie equals `massFlux*Y0`, so a flame near the surface changes the composition on it and not what enters |


and a case needs

```
libs (libburnerBoundaryConditions.so libopticallyThinRadiation.so);
```

in `system/controlDict`. kodes4foam's `tutorials/1D_flat_burner` and
`tutorials/2D_flat_burner` use all of it.

## Name

LM is for low Mach: what the solver is for. It is still the compressible pressure-based
`rhoReactingBuoyantFoam` algorithm underneath, acoustics included.

## License

GPL-3.0, as OpenFOAM; see [LICENSE](LICENSE).
