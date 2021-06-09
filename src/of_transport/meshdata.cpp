#include "meshdata.hpp"

#include <vector>

#include <processorFvPatch.H>

defineTypeNameAndDebug(AuxMeshData, 0);

AuxMeshData::AuxMeshData(const Foam::fvMesh &mesh)
    : Foam::MeshObject<Foam::fvMesh, Foam::GeometricMeshObject, AuxMeshData> (mesh),
      face_to_patch(mesh.nFaces(), -1),
      cell_flags(mesh.nCells(), ECellCategory::regular)
{
    Foam::Info << "Computing face to patch mapping" << Foam::endl;
    for (int patchi = 0; patchi < mesh.boundary().size(); ++patchi) {
        for (int patch_facei = 0; patch_facei < mesh.boundary()[patchi].size(); ++patch_facei) {
            int facei = mesh.boundary()[patchi].start() + patch_facei;
            face_to_patch[facei] = patchi;
        }
    }

    for (int patchi = 0; patchi < mesh.boundary().size(); patchi++) {
        const Foam::fvPatch &p = mesh.boundary()[patchi];
        if (Foam::isA<Foam::processorFvPatch>(p)) {
            for (int celli: p.faceCells()) {
                (uint8_t&)cell_flags[celli] |= (uint8_t)ECellCategory::processor_boundary;
            }
        }
        else {
            for (int celli: p.faceCells()) {
                (uint8_t&)cell_flags[celli] |= (uint8_t)ECellCategory::boundary;
            }
        }
    }
}

AuxMeshData::~AuxMeshData() = default;
