/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    scalarTransportFoam

Group
    grpBasicSolvers

Description
    Passive scalar transport equation solver.

    \heading Solver details
    The equation is given by:

    \f[
        \ddt{T} + \div \left(\vec{U} T\right) - \div \left(D_T \grad T \right)
        = S_{T}
    \f]

    Where:
    \vartable
        T       | Passive scalar
        D_T     | Diffusion coefficient
        S_T     | Source
    \endvartable

    \heading Required fields
    \plaintable
        T       | Passive scalar
        U       | Velocity [m/s]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "meshdata.hpp"
#include "microdomain.hpp"
#include "lsgrad.hpp"

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "surfaceFields.H"

#include <chrono>
#include <vector>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void scalar_transport_explicit_step_internal_faces(
        const volScalarField& Y,
        volScalarField& Y_new,
        const volVectorField& grad_Y,
        const surfaceScalarField& flux,
        const volScalarField& D,
        scalar dt,
        const TMicrodomain& md
        )
{
    const fvMesh& mesh = Y.mesh();
    const auto& owner = mesh.owner();
    const auto& neighbour = mesh.neighbour();
    const auto& V = mesh.V();
    const auto& C = mesh.C();
    const auto& Cf = mesh.Cf();
    const auto& magSf = mesh.magSf();
    const auto& weights = mesh.weights();
    const auto& NODC = mesh.nonOrthDeltaCoeffs();
    const auto& NOCV = mesh.nonOrthCorrectionVectors();

    for (label celli = md.start_cell; celli < md.end_cell; celli++) {
        Y_new[celli] = Y[celli];
    }

    for (label facei = md.start_face; facei < md.end_face; facei++) {
        label own = owner[facei];
        label nei = neighbour[facei];
        // convective term
        label upwind_cell = flux[facei] > 0.0 ? own : nei;
        vector cell_to_face = Cf[facei] - C[upwind_cell];
        scalar Y_f = Y[upwind_cell] + (cell_to_face & grad_Y[upwind_cell]);
        // diffusive term
        scalar w  = weights[facei];
        scalar interp_D          = w * D[own]      + (1.0 - w) * D[nei];
        vector interp_grad_Y     = w * grad_Y[own] + (1.0 - w) * grad_Y[nei];
        // equivalent to "corrected"
        scalar snGrad_Y = (Y[nei] - Y[own]) * NODC[facei]
                + (NOCV[facei] & interp_grad_Y);

        scalar f = magSf[facei] * (Y_f * flux[facei] - snGrad_Y * interp_D);
        Y_new[own] -= f * dt / V[own];
        Y_new[nei] += f * dt / V[nei];
    }
}

void scalar_transport_explicit_step_boundary_faces(
        const volScalarField& Y,
        volScalarField& Y_new,
        const volVectorField& grad_Y,
        const surfaceScalarField& flux,
        const volScalarField& D,
        scalar dt
        )
{
    const fvMesh& mesh = Y.mesh();
    const auto& V = mesh.V();
    const auto& C = mesh.C();
    const auto& weights = mesh.weights();
    const auto& NODC = mesh.nonOrthDeltaCoeffs();
    const auto& NOCV = mesh.nonOrthCorrectionVectors();

    for (int patchi = 0; patchi < mesh.boundary().size(); ++patchi) {
        const fvPatch &p = mesh.boundary()[patchi];
        const auto& p_face_cells = p.faceCells();
        const auto& p_weights = weights.boundaryField()[patchi];
        const auto& pFlux = flux.boundaryField()[patchi];
        if (p.coupled()) {
            const scalarField pYRight(Y.boundaryField()[patchi].patchNeighbourField());
            const vectorField pGradYRight(grad_Y.boundaryField()[patchi].patchNeighbourField());
            const vectorField pCRight(C.boundaryField()[patchi].patchNeighbourField());
            const scalarField pDRight(D.boundaryField()[patchi].patchNeighbourField());

            const auto& h = NODC.boundaryField()[patchi];
            const auto& hv= NOCV.boundaryField()[patchi];

            for (int patch_facei = 0; patch_facei < p.size(); ++patch_facei) {
                label celli = p_face_cells[patch_facei];
                // convective flux:
                scalar flux = pFlux[patch_facei];
                scalar Y_f = 0.0;
                if (flux > 0.0) {
                    // outflow
                    vector cell_to_face = p.Cf()[patch_facei] - C[celli];
                    Y_f = Y[celli] + (cell_to_face & grad_Y[celli]);
                }
                else {
                    //inflow
                    vector cell_to_face = p.Cf()[patch_facei] - pCRight[patch_facei];
                    Y_f = pYRight[patch_facei] + (cell_to_face & pGradYRight[patch_facei]);
                }
                // diffusive flux:
                scalar w = p_weights[patch_facei];
//                scalar interp_rho    = w * prhoLeft[patch_facei]   + (1.0 - w) * prhoRight[patch_facei];
                scalar interp_D      = w * D[celli]     + (1.0 - w) * pDRight[patch_facei];
                vector interp_grad_Y = w * grad_Y[celli] + (1.0 - w) * pGradYRight[patch_facei];
                scalar snGrad_Y = (pYRight[patch_facei] - Y[celli]) * h[patch_facei] + (hv[patch_facei] & interp_grad_Y);
                scalar f = (Y_f * flux - snGrad_Y * interp_D) * p.magSf()[patch_facei];
                Y_new[celli] -= f * dt / V[celli];
            }
        }
        else {
            scalarField p_snGrad_Y(Y.boundaryField()[patchi].snGrad());
            const auto& p_D = D.boundaryField()[patchi];
            const auto& p_Y = Y.boundaryField()[patchi];
            for (int patch_facei = 0; patch_facei < p.size(); ++patch_facei) {
                label celli = p_face_cells[patch_facei];
                scalar flux = pFlux[patch_facei];
                scalar Y_f = 0.0;
                if (flux > 0.0) {
                    // outflow
                    vector cell_to_face = p.Cf()[patch_facei] - C[celli];
                    Y_f = Y[celli] + (cell_to_face & grad_Y[celli]);
                }
                else {
                    //inflow
                    Y_f = p_Y[patch_facei];
                }
                scalar interp_D = p_D[patch_facei];
                scalar snGrad_Y = p_snGrad_Y[patch_facei];
                scalar f = (Y_f * flux - snGrad_Y * interp_D) * p.magSf()[patch_facei];
                Y_new[celli] -= f * dt / V[celli];
            }
        }
    }
}

static inline double limit_f_venkatakrishnan(double delta_max, double delta_min, double extrapolate, double eps_sqr){
    double l = 1.0;
    if (extrapolate > SMALL) {
        l = (delta_max * delta_max + eps_sqr + 2.0 * delta_max * extrapolate )/
                (delta_max * delta_max + 2.0 * extrapolate * extrapolate + extrapolate * delta_max + eps_sqr);
    }
    else if (extrapolate < -SMALL) {
        l = (delta_min * delta_min + eps_sqr + 2.0 * delta_min * extrapolate )/
                (delta_min * delta_min + 2.0 * extrapolate * extrapolate + extrapolate * delta_min + eps_sqr);
    }
    return l;
}

static void apply_scalar_limiter(volVectorField &grad,
                          const volScalarField &val,
                          double limiter_eps_sqr = 1e-04)
{
    const Foam::fvMesh &mesh = grad.mesh();
    const auto& owner = mesh.owner();
    const auto& neighbour = mesh.neighbour();
    const auto& C = mesh.C();
    const auto& points = mesh.points();
    const auto& meshdata = AuxMeshData::New(mesh);
    for(int celli = 0 ; celli < mesh.nCells() ; celli ++){
        const Foam::cell &c = mesh.cells()[celli];
        double val_own = val[celli];
        double val_min = val_own;
        double val_max = val_own;
        for(int j = 0; j < c.nFaces(); ++j) {
            int facei = c[j];
            if (mesh.isInternalFace(facei)) {
                label other = owner[facei];
                if (other == celli) {
                    other = neighbour[facei];
                }
                double val_other = val[other];
                if(val_other < val_min){
                    val_min = val_other;
                }
                if(val_other > val_max){
                    val_max = val_other;
                }
            }
            else{
                label patchi  = meshdata.face_to_patch[facei];
                if (patchi >= 0){
                    label   patch_face_i = mesh.boundaryMesh()[patchi].whichFace(facei);
                    double  val_other    = val.boundaryField()[patchi][patch_face_i];
                    if(val_other < val_min){
                        val_min = val_other;
                    }
                    if(val_other > val_max){
                        val_max = val_other;
                    }
                }
            }
        }

        double delta_min = val_min - val_own;
        double delta_max = val_max - val_own;
        double eps_sqr = limiter_eps_sqr * (val_min * val_min + val_max * val_max);
        vector grad_own = grad[celli];
        double l = 1.0;

        const labelList &plist = mesh.cellPoints()[celli];
        for (int j = 0; j < plist.size(); ++j) {
            vector cell_to_point = points[plist[j]] - C[celli];
            double extrapolate  = (cell_to_point & grad_own);
            l = min(l , limit_f_venkatakrishnan(delta_max, delta_min, extrapolate, eps_sqr));
        }
        grad[celli] *= l;
    }
}

static void apply_scalar_limiter_minmod_scalar(
        volVectorField &grad,
        const volScalarField &val,
        const AuxMeshData& mesh_data,
        const TMicrodomain& md)
{
    const fvMesh &mesh = grad.mesh();

    for (label celli = md.start_cell; celli < md.end_cell; ++celli) {
        const cell &c = mesh.cells()[celli];
        scalar lim = 1.0;
        for (label facei: c) {
            if (mesh.isInternalFace(facei)) {
                label other = mesh.owner()[facei];
                if (other == celli) {
                    other = mesh.neighbour()[facei];
                }

                vector cc = mesh.C()[other] - mesh.C()[celli];
                scalar extrapolate = cc & grad[celli];

                scalar val_diff = val[other] - val[celli];
                if (extrapolate != 0.0) {
                    lim = min(lim, max(0.0, val_diff / extrapolate));
                }
            }
            else {
                label patchi = mesh_data.face_to_patch[facei];
                if (patchi < 0){
                    continue;
                }
                const fvPatch &patch = mesh.boundary()[patchi];
                label patch_face_i = patch.patch().whichFace(facei);

                vector cc = patch.Cf()[patch_face_i] - mesh.C()[celli];
                scalar extrapolate = cc & grad[celli];

                scalar val_diff = val.boundaryField()[patchi][patch_face_i] - val[celli];
                if (extrapolate != 0.0) {
                    lim = min(lim, max(0.0, val_diff / extrapolate));
                }
            }
        }

        grad[celli] = grad[celli] * lim;
    }
}

inline void limit_face
(
    scalar& limiter,
    const scalar maxDelta,
    const scalar minDelta,
    const scalar extrapolate
)
{
    scalar r = 1;

    if (extrapolate > SMALL)
    {
        r = maxDelta/extrapolate;
    }
    else if (extrapolate < -SMALL)
    {
        r = minDelta/extrapolate;
    }
    else
    {
        return;
    }

    limiter = min(limiter, r);
}

static void apply_scalar_limiter_minmod_2 (
        volVectorField &grad,
        const volScalarField &val,
//        std::vector<double> max,
//        std::vector<double> min,
//        std::vector<double> lim,
        const AuxMeshData& mesh_data,
        const TMicrodomain& md)
{
    const fvMesh& mesh = grad.mesh();

    const auto& C = mesh.C();
    const auto& Cf = mesh.Cf();

//    for (int facei = md.start_face; facei < md.end_face; facei++) {
//        int own = mesh.owner()[facei];
//        int nei = mesh.neighbour()[facei];

//        if (val[nei] > max[own]) {
//            max[own] = val[nei];
//        }
//        if (val[nei] < min[own]) {
//            min[own] = val[nei];
//        }
//        if (val[own] > max[nei]) {
//            max[nei] = val[own];
//        }
//        if (val[own] < min[nei]) {
//            min[nei] = val[own];
//        }
//    }

    for (int celli = md.start_cell; celli < md.end_cell; ++celli) {
        scalar v = val[celli];
        scalar min_val = v;
        scalar max_val = v;

        for (int nei : mesh.cellCells()[celli]) {
            scalar v_nei = val[nei];
            if (v_nei > max_val) {
                max_val = v_nei;
            }
            if (v_nei < min_val) {
                min_val = v_nei;
            }
        }

        if (mesh_data.cell_flags[celli] & (ECellCategory::boundary | ECellCategory::processor_boundary)) {
            for (label facei: mesh.cells()[celli]) {
                if (mesh.isInternalFace(facei)) {
                    continue;
                }
                label patchi = mesh_data.face_to_patch[facei];
                label patch_face_i = mesh.boundaryMesh()[patchi].whichFace(facei);
                scalar v_nei = val.boundaryField()[patchi][patch_face_i];
                if (v_nei > max_val) {
                    max_val = v_nei;
                }
                if (v_nei < min_val) {
                    min_val = v_nei;
                }
            }
        }

        max_val -= v;
        min_val -= v;

        vector c = C[celli];
        scalar limiter = 1.0;
        for (int nei : mesh.cellCells()[celli]) {
            vector cc = mesh.C()[nei] - c;
            scalar extrapolate = cc & grad[celli];
            if (limiter * extrapolate > max_val) {
                limiter = max_val / extrapolate;
            }
            if (limiter * extrapolate < min_val) {
                limiter = min_val / extrapolate;
            }
        }

        if (mesh_data.cell_flags[celli] & (ECellCategory::boundary | ECellCategory::processor_boundary)) {
            for (label facei: mesh.cells()[celli]) {
                if (mesh.isInternalFace(facei)) {
                    continue;
                }
                label patchi = mesh_data.face_to_patch[facei];
                const fvPatch& p = mesh.boundary()[patchi];
                label patch_face_i = p.patch().whichFace(facei);
                vector cc = p.Cf()[patch_face_i] - c;
                scalar extrapolate = cc & grad[celli];
                if (limiter * extrapolate > max_val) {
                    limiter = max_val / extrapolate;
                }
                if (limiter * extrapolate < min_val) {
                    limiter = min_val / extrapolate;
                }
            }
        }

        grad[celli] *= limiter;
    }
}

int main(int argc, char *argv[])
{
    argList::addOption("method", "M", "Method to use for species transport: OpenFOAM, singlePass, microdomains, microdomainsOpenMP");
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    Info << "Creating auxiliary mesh data" << endl;

    const AuxMeshData& mesh_data = AuxMeshData::New(mesh);

    word method = args.get<word>("method");

    std::vector<TMicrodomain> microdomains;
    std::vector<label> microdomain;
    if (method.starts_with("microdomains")) {
        Info << "Reading microdomains" << endl;
        read_microdomains(microdomain, microdomains, runTime, mesh);
    }
    TMicrodomain whole_mesh{.start_face=0, .end_face=mesh.nInternalFaces(), .start_cell=0, .end_cell=mesh.nCells()};

    Info << "Creating least squares grad scheme" << Foam::endl;

    TLeastSquaresGrad lsgrad(mesh);

    Info<< "Creating specie fields" << endl;

    const size_t n_species = 10;
    std::vector<volScalarField> Y;
    std::vector<volScalarField> Ynew;
    std::vector<volVectorField> gradY;
    Y.reserve(n_species);
    Ynew.reserve(n_species);
    gradY.reserve(n_species);
    for (size_t s = 0; s < n_species; ++s) {
        std::string name = std::string("Y_") + std::to_string(s);
        Info<< "Creating specie field " << name << endl;
        Y.emplace_back(IOobject(name, runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                             mesh, dimensionedScalar(name, dimless, 0.0));
        Ynew.emplace_back(IOobject(name+"new", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                             mesh, dimensionedScalar(name, dimless, 0.0));
        gradY.emplace_back(IOobject("grad"+name, runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                             mesh, dimensionedVector(name, dimless/dimLength, Zero));
    }

    volScalarField D    (IOobject("D",     runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), mesh,
                         dimensionedScalar("", dimArea/dimTime, 1.0));

//    std::vector<label> md(mesh.nCells());
//    for (label celli = 0; celli < mesh.nCells(); celli++) {
//        md[celli] = microdomain[celli];
//    }

//    Info << "Sorted: " << std::is_sorted(md.begin(), md.end()) << endl;

    Info<< "Reading field U\n" << endl;
    volVectorField U(IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh);

    Info<< "Reading/calculating face flux field phi\n" << endl;
    surfaceScalarField phi(IOobject("phi", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
                           fvc::flux(U)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        auto start = std::chrono::steady_clock::now();

        Info << "Solving species " << endl;
        if (method == "OpenFOAM") {
            // Option 1: standard OpenFOAM equations
            for (size_t s = 0; s < Y.size(); ++s) {
                solve
                        (
                            fvm::ddt(Y[s])
                            ==
                            - fvc::div(phi, Y[s])
                            + fvc::laplacian(D, Y[s])
                            );
            }
        }
        else if (method == "singlePass") {
            for (size_t s = 0; s < Y.size(); ++s) {
                gradY[s] = fvc::grad(Y[s]);
                apply_scalar_limiter_minmod_2(gradY[s], Y[s], mesh_data, whole_mesh);
                scalar_transport_explicit_step_internal_faces(Y[s], Ynew[s], gradY[s], phi, D, runTime.deltaTValue(), whole_mesh);
                scalar_transport_explicit_step_boundary_faces(Y[s], Ynew[s], gradY[s], phi, D, runTime.deltaTValue());
                Ynew[s].correctBoundaryConditions();
                Y[s] = Ynew[s];
            }
        }
        else if (method == "microdomains") {
            for (size_t s = 0; s < Y.size(); ++s) {
                gradY[s] = fvc::grad(Y[s]);
                apply_scalar_limiter_minmod_2(gradY[s], Y[s], mesh_data, whole_mesh);
                for (const auto& md : microdomains) {
                    scalar_transport_explicit_step_internal_faces(Y[s], Ynew[s], gradY[s], phi, D, runTime.deltaTValue(), md);
                }
                scalar_transport_explicit_step_boundary_faces(Y[s], Ynew[s], gradY[s], phi, D, runTime.deltaTValue());
                Ynew[s].correctBoundaryConditions();
                Y[s] = Ynew[s];
            }
        }
        else if (method == "microdomains_lsgrad") {
            for (size_t s = 0; s < Y.size(); ++s) {
                lsgrad.grad(Y[s], gradY[s], whole_mesh);
                apply_scalar_limiter_minmod_2(gradY[s], Y[s], mesh_data, whole_mesh);
                for (const auto& md : microdomains) {
                    scalar_transport_explicit_step_internal_faces(Y[s], Ynew[s], gradY[s], phi, D, runTime.deltaTValue(), md);
                }
                scalar_transport_explicit_step_boundary_faces(Y[s], Ynew[s], gradY[s], phi, D, runTime.deltaTValue());
                Ynew[s].correctBoundaryConditions();
                Y[s] = Ynew[s];
            }
        }
        else if (method == "microdomains_lsgrad_2") {
            for (size_t s = 0; s < Y.size(); ++s) {
                for (const auto& md : microdomains) {
                    lsgrad.grad(Y[s], gradY[s], md);
                    apply_scalar_limiter_minmod_2(gradY[s], Y[s], mesh_data, md);
                    scalar_transport_explicit_step_internal_faces(Y[s], Ynew[s], gradY[s], phi, D, runTime.deltaTValue(), md);
                }
                scalar_transport_explicit_step_boundary_faces(Y[s], Ynew[s], gradY[s], phi, D, runTime.deltaTValue());
                Ynew[s].correctBoundaryConditions();
                Y[s] = Ynew[s];
            }
        }
        else if (method == "microdomains_lsgrad_swapped") {
            for (const auto& md : microdomains) {
                for (size_t s = 0; s < Y.size(); ++s) {
                    lsgrad.grad(Y[s], gradY[s], md);
                    apply_scalar_limiter_minmod_2(gradY[s], Y[s], mesh_data, md);
                    scalar_transport_explicit_step_internal_faces(Y[s], Ynew[s], gradY[s], phi, D, runTime.deltaTValue(),
                                                                  md);
                }
            }
            for (size_t s = 0; s < Y.size(); ++s) {
                scalar_transport_explicit_step_boundary_faces(Y[s], Ynew[s], gradY[s], phi, D, runTime.deltaTValue());
                Ynew[s].correctBoundaryConditions();
                Y[s] = Ynew[s];
            }
        }
        else {
            FatalError << "Unknown method: " << method << abort(FatalError);
        }

        Info << "Elapsed clock time = " << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() << endl;
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}
