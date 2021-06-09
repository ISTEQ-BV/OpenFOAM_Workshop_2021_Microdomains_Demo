#include "lsgrad.hpp"

#include "fvc.H"
#include "argList.H"

using namespace Foam;

template <typename T, typename Expected>
void print_cell_error(const Foam::word &name, const GeometricField<T, fvPatchField, volMesh> &f, const Expected &expected) {
    const GeometricField<T, fvPatchField, volMesh> error(f - expected);

    scalar error_max = 0.0;
    scalar error_sum_sqr = 0.0;
    for (int celli = 0; celli < f.mesh().nCells(); celli++) {
        if (mag(error[celli]) > error_max) {
            error_max = mag(error[celli]);
        }
        error_sum_sqr += magSqr(error[celli]);
    }
    Foam::reduce(error_max, Foam::maxOp<scalar>());
    Foam::reduce(error_sum_sqr, Foam::sumOp<scalar>());
    Foam::Info << name << '\t' << error_sum_sqr << '\t' << error_max << Foam::endl;
}

void test_lsgrad_zero(const fvMesh &mesh) {

    Info << "Computing weights" << endl;

    TLeastSquaresGrad lsgrad(mesh);

    volScalarField fs(IOobject("fs", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                      mesh, dimensionedScalar("", dimless, 0.0), "zeroGradient");

    volVectorField fs_grad(IOobject("grad(fs)", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                               mesh, dimensionedVector("", dimless/dimLength, Zero));

    Info << "Computing grad" << endl;
    TMicrodomain whole_mesh{.start_face=0, .end_face=mesh.nInternalFaces(), .start_cell=0, .end_cell=mesh.nCells()};
    lsgrad.grad(fs, fs_grad, whole_mesh);

    print_cell_error("zero", fs_grad, dimensionedVector("", dimless/dimLength, Zero));
}

void test_lsgrad_const(const fvMesh &mesh) {

    Info << "Computing weights" << endl;

    TLeastSquaresGrad lsgrad(mesh);

    volScalarField fs(IOobject("fs", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                      mesh, dimensionedScalar("", dimless, 1.0), "zeroGradient");

    volVectorField fs_grad(IOobject("grad(fs)", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                               mesh, dimensionedVector("", dimless/dimLength, Zero));

    Info << "Computing grad" << endl;
    TMicrodomain whole_mesh{.start_face=0, .end_face=mesh.nInternalFaces(), .start_cell=0, .end_cell=mesh.nCells()};
    lsgrad.grad(fs, fs_grad, whole_mesh);

    print_cell_error("const", fs_grad, dimensionedVector("", dimless/dimLength, Zero));
}

void test_lsgrad_linear(const fvMesh &mesh) {

    Info << "Computing weights" << endl;

    TLeastSquaresGrad lsgrad(mesh);

    volScalarField fs(IOobject("fs", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                      mesh.C().component(1));

    volVectorField fs_grad(IOobject("grad(fs)", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                               mesh, dimensionedVector("", dimless/dimLength, Zero));

    Info << "Computing grad" << endl;

    TMicrodomain whole_mesh{.start_face=0, .end_face=mesh.nInternalFaces(), .start_cell=0, .end_cell=mesh.nCells()};
    lsgrad.grad(fs, fs_grad, whole_mesh);

    print_cell_error("linear", fs_grad, dimensionedVector("", dimless/dimLength, {0.0, 1.0, 0.0}));
}

int main(int argc, char *argv[]) {
    /* Set root case */
    Foam::argList args(argc, argv);
    if (!args.checkRootCase()) {
        Foam::FatalError.exit();
    }
    /* Create time */
    Foam::Info<< "Create time\n" << Foam::endl;
    Foam::Time runTime(Foam::Time::controlDictName, args);
    /* Create mesh */
    Foam::Info  << "Create mesh for time = "  << runTime.timeName() << Foam::nl << Foam::endl;
    Foam::fvMesh mesh ( Foam::IOobject ( Foam::fvMesh::defaultRegion, runTime.timeName(), runTime, Foam::IOobject::MUST_READ ) );

    test_lsgrad_zero(mesh);
    test_lsgrad_const(mesh);
    test_lsgrad_linear(mesh);

    return 0;
}
