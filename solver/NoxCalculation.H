basicMultiComponentMixture& composition = thermo.composition();


volScalarField cNOx =   YNOx*rho/NOe.W()/1000;
volScalarField cOH  =   composition.Y("OH")*rho/OHe.W()/1000;
volScalarField cN2  =   composition.Y("N2")*rho/N2e.W()/1000;
volScalarField cO2  =   composition.Y("O2")*rho/O2e.W()/1000;
volScalarField cO   =   composition.Y("O")*rho/Oe.W()/1000;
volScalarField cH   =   composition.Y("H")*rho/He.W()/1000;



forAll(YNOx,celli)
{
scalar k1pos = 7.6e13 * ::exp(-38000/T[celli]);
scalar k2pos = 6.4e9 * T[celli] * ::exp(-3150/T[celli]);
scalar k3pos = 4.1e13;
scalar k1neg = 1.6e13;
scalar k2neg = 1.5e9 * T[celli] * ::exp(-19500/T[celli]);
scalar K = (k1pos/k1neg)*(k2pos/k2neg);
scalar A = 2.0 * k1pos*cO[celli]*cN2[celli];
scalar B = 1.0 - ::pow(cNOx[celli],2.0)/(K*cO2[celli]*cN2[celli]);
scalar C = 1.0 + k1neg*cNOx[celli]/(k2pos*cO2[celli]+k3pos*cOH[celli]);
RRNOx[celli] = A*B/C;
}

tmp<fv::convectionScheme<scalar> > mvConvection
(
fv::convectionScheme<scalar>::New
(
mesh,
fields,
phi,
mesh.divScheme("div(phi,Yi_h)")
)
);


solve
(
fvm::ddt(rho, YNOx)
+ mvConvection->fvmDiv(phi, YNOx)
- fvm::laplacian(thermo.mu() + rSct * turbulence->mut(), YNOx)
==
RRNOx,
mesh.solver("Yi")
);
