#include "microdomain.hpp"

#include "volFields.H"

using Foam::label;

void read_microdomains(std::vector<Foam::label> &microdomain, std::vector<TMicrodomain>& microdomains, const Foam::Time &runTime, const Foam::fvMesh &mesh) {
    Foam::volScalarField microdomain_field(Foam::IOobject("microdomain", runTime.timeName(), mesh, Foam::IOobject::MUST_READ, Foam::IOobject::NO_WRITE), mesh);

    microdomain.assign(mesh.nCells(), -1);

    for (label celli = 0; celli < mesh.nCells(); celli++) {
        label m = static_cast<label>(microdomain_field[celli]);
        microdomain[celli] = m;
        if (m+1 == microdomains.size()) {
            continue;
        }
        else if (m == microdomains.size()) {
//            Info << "New microdomain " << m << endl;
            if (m > 0)
                microdomains.back().end_cell = celli;
            microdomains.push_back({-1, -1, celli, -1});
        }
        else {
            Foam::FatalError << "Microdomains are not ordered properly" << Foam::abort(Foam::FatalError);
        }
    }
    microdomains.back().end_cell = mesh.nCells();

    std::vector<std::vector<label>> microdomain_faces;
    microdomain_faces.resize(microdomains.size());

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++) {
        label own = mesh.owner()[facei];
        label md = microdomain[own];
        microdomain_faces[md].push_back(facei);
    }

    for (size_t i = 0; i < microdomains.size(); ++i) {
        auto & md = microdomains[i];
        const auto & faces = microdomain_faces[i];
        label n_faces = faces.size();
        auto p = std::minmax_element(faces.begin(), faces.end());
        label start_face = *p.first;
        label end_face = *p.second + 1;
        if (end_face - start_face != n_faces) {
            Foam::FatalError << "Faces for microdomains are not contiguous" << Foam::abort(Foam::FatalError);
        }
        md.start_face = start_face;
        md.end_face = end_face;
    }
};
