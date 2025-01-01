volScalarField YNOx
(
IOobject
(
"YNOx",
runTime.timeName(),
mesh,
IOobject::READ_IF_PRESENT,
IOobject::AUTO_WRITE
),
mesh,
dimensionedScalar("zero", dimless, 0.0),
zeroGradientFvPatchScalarField::typeName
);
volScalarField RRNOx
(
IOobject
(
"RRNOx",
runTime.timeName(),
mesh,
IOobject::NO_READ,
IOobject::NO_WRITE
),
mesh,
dimensionedScalar("zero", dimensionSet (1,-3,-1,0,0), 0.0),
zeroGradientFvPatchScalarField::typeName
);
fields.add(YNOx);
fields.add(RRNOx);



dictionary NOXDATA(thermo.subDict("NOxData"));

gasHThermoPhysics NOe(NOXDATA.subDict("NO"));
gasHThermoPhysics OHe(NOXDATA.subDict("OH"));
gasHThermoPhysics N2e(NOXDATA.subDict("N2"));
gasHThermoPhysics O2e(NOXDATA.subDict("O2"));
gasHThermoPhysics Ne(NOXDATA.subDict("N"));
gasHThermoPhysics Oe(NOXDATA.subDict("O"));
gasHThermoPhysics He(NOXDATA.subDict("H"));

dimensionedScalar Sct
(
combustionProperties.lookup("Sct")
);
dimensionedScalar rSct = 1.0 / Sct;
