#ifndef MESHDATA_HPP
#define MESHDATA_HPP

#include "MeshObject.H"
#include "fvMesh.H"

#include <vector>

using Foam::word;

enum ECellCategory : uint8_t {
    /// Regular category, not adjacent to processor boundaries
    regular                = 0x0,
    /// Cell is directly adjacent to propcessor boundary
    processor_boundary     = 0x1,
    /// A cell that near the nonprocessor boundary
    boundary               = 0x2,
};

class AuxMeshData
    : public Foam::MeshObject<Foam::fvMesh, Foam::GeometricMeshObject, AuxMeshData>
{
public:
    TypeName("MeshData")
    explicit AuxMeshData(const Foam::fvMesh &mesh);

    virtual ~AuxMeshData();

public:
    std::vector<int> face_to_patch;
    std::vector<ECellCategory> cell_flags;
};

#endif // MESHDATA_HPP
