Info<< "Reading field T\n" << endl;

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
T.oldTime();

volVectorField U
(
    IOobject
    (
        "U",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("U", dimVelocity, vector::zero),
    "zeroGradient"
);

surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    linearInterpolate(U)
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    Uf & mesh.Sf()
);

// Courant numbers to write out
volScalarField Co
(
    IOobject("Co", runTime.constant(), mesh),
    CourantNo(phi, dt)
);
Co.write();
// Make a indicator function based on the CFL number of a cell
volScalarField INDI("INDI", T);
volScalarField INDE("INDE", T);
forAll(mesh.C(), celli)
{
    if (mag(Co[celli]) > CoLimit)
    {
        INDI[celli] = 1;
        INDE[celli] = 0;
    }
    else
    {
        INDI[celli] = 0;
        INDE[celli] = 1;
    }
};
INDE.write();
// Intermediate values of T
volScalarField KK("KK", T);

