#include "lsgrad.hpp"
#include "meshdata.hpp"
#include "microdomain.hpp"

#include <cmath>
#include <vector>

#include <vector.H>
#include <fvMesh.H>
#include <volFields.H>
#include <processorFvPatch.H>

static void compute_ls_coeffs_for_cell(
        int celli,
        std::vector<Foam::label>&  cell_idx,
        std::vector<Foam::vector>& cell_coeff,
        std::vector<Foam::label>&  face_idx,
        std::vector<Foam::label>&  face_patch,
        std::vector<Foam::vector>& face_coeff,
        const Foam::fvMesh &mesh,
        const AuxMeshData& mesh_data)
{
    Foam::symmTensor R = Foam::Zero;

    Foam::vector c = mesh.C()[celli];

    for (int nei: mesh.cellCells()[celli]) {
        Foam::vector r = mesh.C()[nei] - c;
        double w2 = 1.0 / Foam::magSqr(r);
        R += w2 * sqr(r);
    }

    if ((uint8_t)mesh_data.cell_flags[celli] & ((uint8_t)ECellCategory::boundary | (uint8_t)ECellCategory::processor_boundary)) {
        for (Foam::label facei: mesh.cells()[celli]) {
            if (! mesh.isInternalFace(facei)) {
                int patchi = mesh_data.face_to_patch[facei];
                const Foam::polyPatch& pp = mesh.boundaryMesh()[patchi];
                Foam::label patch_face_i = pp.whichFace(facei);
                Foam::vector r = mesh.C().boundaryField()[patchi][patch_face_i] - c;
                double w2 = 1.0 / Foam::magSqr(r);
                R += w2 * sqr(r);
            }
        }
    }

    Foam::symmTensor invR = inv(R);

    for (int nei: mesh.cellCells()[celli]) {
            Foam::vector r = mesh.C()[nei] - c;
            double w2 = 1.0 / Foam::magSqr(r);
            Foam::vector coeff = w2 * (invR & r);
            cell_idx.push_back(nei);
            cell_coeff.push_back(coeff);
    }

    if ((uint8_t)mesh_data.cell_flags[celli] & ((uint8_t)ECellCategory::boundary | (uint8_t)ECellCategory::processor_boundary)) {
        for (Foam::label facei: mesh.cells()[celli]) {
            if (! mesh.isInternalFace(facei)) {
                int patchi = mesh_data.face_to_patch[facei];
                const Foam::polyPatch& pp = mesh.boundaryMesh()[patchi];
                Foam::label patch_face_i = pp.whichFace(facei);
                Foam::vector r = mesh.C().boundaryField()[patchi][patch_face_i] - c;
                double w2 = 1.0 / Foam::magSqr(r);
                Foam::vector coeff = w2 * (invR & r);
                face_idx.push_back(patch_face_i);
                face_patch.push_back(patchi);
                face_coeff.push_back(coeff);
            }
        }
    }
}

TLeastSquaresGrad::TLeastSquaresGrad(const Foam::fvMesh &mesh)
{
    compute_coeffs(mesh);
}

void TLeastSquaresGrad::compute_coeffs(const Foam::fvMesh &mesh)
{
    const auto& mesh_data = AuxMeshData::New(mesh);

    cell_row.reserve(mesh.nCells()+1);
    face_row.reserve(mesh.nCells()+1);

    cell_idx  .reserve(2*mesh.nInternalFaces());
    cell_coeff.reserve(2*mesh.nInternalFaces());
    face_patch.reserve(mesh.nBoundaryFaces());
    face_idx  .reserve(mesh.nBoundaryFaces());
    face_coeff.reserve(mesh.nBoundaryFaces());

    cell_row.push_back(0);
    face_row.push_back(0);
    for (Foam::label celli = 0; celli < mesh.nCells(); celli++) {
        compute_ls_coeffs_for_cell(celli, cell_idx, cell_coeff, face_idx, face_patch, face_coeff, mesh, mesh_data);
        cell_row.push_back(cell_idx.size());
        face_row.push_back(face_idx.size());
    }
}

void TLeastSquaresGrad::grad(const Foam::volScalarField &src, Foam::volVectorField &dst, const TMicrodomain &md)
{
    for (int celli = md.start_cell; celli < md.end_cell; celli++) {
        Foam::vector grad = Foam::Zero;

        int32_t cell_start = cell_row[celli  ];
        int32_t cell_end   = cell_row[celli+1];
        for(int i = cell_start; i < cell_end; i++){
            double delta = (src[cell_idx[i]] - src[celli]);
            grad += cell_coeff[i] * delta;
        }

        int32_t face_start = face_row[celli  ];
        int32_t face_end   = face_row[celli+1];
        for(int i = face_start; i < face_end; i++) {
            int pfacei = face_idx[i];
            int patchi = face_patch[i];
            const auto& p = src.boundaryField()[patchi];
            double delta = (p[pfacei] - src[celli]);
            grad += face_coeff[i] * delta;
        }

        dst[celli] = grad;
    }
}
