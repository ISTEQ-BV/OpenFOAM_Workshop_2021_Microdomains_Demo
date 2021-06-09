#ifndef LSGRAD_HPP
#define LSGRAD_HPP

#include "microdomain.hpp"

#include "fvMesh.H"

#include <vector>
#include <cstdint>

struct TLeastSquaresGrad
{
    TLeastSquaresGrad(const Foam::fvMesh& mesh);
    void compute_coeffs(const Foam::fvMesh& mesh);
    void grad(const Foam::volScalarField& src, Foam::volVectorField& dst, const TMicrodomain& md);

    std::vector<int32_t>  cell_row;
    std::vector<int32_t>  cell_idx;
    std::vector<Foam::vector>   cell_coeff;
    std::vector<int32_t>  face_row;
    std::vector<int32_t>  face_patch;
    std::vector<int32_t>  face_idx;
    std::vector<Foam::vector>   face_coeff;
};

#endif // LSGRAD_HPP
