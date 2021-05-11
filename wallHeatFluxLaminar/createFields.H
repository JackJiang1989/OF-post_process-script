Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
IOobject
(
    "transportProperties",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
)
);

dimensionedScalar tc
(
"tc",
transportProperties.lookup("tc")
);


volScalarField T
(
IOobject
(
    "T",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
),
mesh
    );
