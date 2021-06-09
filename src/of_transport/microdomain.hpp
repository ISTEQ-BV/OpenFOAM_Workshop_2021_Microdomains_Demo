#ifndef MICRODOMAIN_HPP
#define MICRODOMAIN_HPP

#include "Time.H"
#include "fvMesh.H"

#include "vector"

struct TMicrodomain {
    int start_face;
    int end_face;
    int start_cell;
    int end_cell;
};

void read_microdomains(std::vector<Foam::label>& microdomain, std::vector<TMicrodomain>& microdomains, const Foam::Time& runTime, const Foam::fvMesh& mesh);

#endif // MICRODOMAIN_HPP
