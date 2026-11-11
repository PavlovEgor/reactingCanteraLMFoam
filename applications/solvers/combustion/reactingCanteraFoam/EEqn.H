{
    volScalarField& he = thermo.he();

    radiation->correct();

    // Heat conduction is - laplacian(kappa/Cp, hs) less the part of it that
    // is a gradient of composition; that part and the enthalpy the diffusive
    // fluxes carry are in hsFlux
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      - dpdt
      - fvm::laplacian(transport.alpha(), he, "laplacian(alpha,h)")
      + fvc::div(transport.hsFlux())
     ==
        Qdot
      + radiation->Sh(thermo, he)
      + fvOptions(rho, he)
    );

    if (MRF.active())
    {
        EEqn += fvc::div(MRF.phi(), p);
    }

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
