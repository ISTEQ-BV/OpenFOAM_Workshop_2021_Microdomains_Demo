#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"

#include <algorithm>
#include <vector>
#include <numeric>
#include <fstream>

void save_list(const std::string &filename, const std::vector<int> &cell_order)
{
    std::ofstream cell_order_file(filename);
    cell_order_file << R"(/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v1712                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       labelList;
    location    "constant";
    object      xxx;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

(
)";
    for (int celli : cell_order) {
        cell_order_file << celli << '\n';
    }
    cell_order_file << ");\n";
}

int main(int argc, char **argv) {
    Foam::argList::addOption("nodes", "N", "number of MPI nodes");
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

    Foam::volScalarField cellDist   (Foam::IOobject("cellDist",    runTime.timeName(), mesh, Foam::IOobject::MUST_READ), mesh);
    Foam::volScalarField node       (Foam::IOobject("node",        runTime.timeName(), mesh, Foam::IOobject::NO_READ),   mesh, Foam::dimless);
    Foam::volScalarField microdomain(Foam::IOobject("microdomain", runTime.timeName(), mesh, Foam::IOobject::NO_READ),   mesh, Foam::dimless);

    // number of MPI nodes
    int n_nodes = args.getOrDefault("nodes", 1);
    Foam::Info << "Number of nodes = " << n_nodes << Foam::endl;

    // Total number of microdomains
    int n_domains = *std::max_element(cellDist.begin(), cellDist.end()) + 1;
    Foam::Info << "Number of domains = " << n_domains << Foam::endl;

    std::vector<int> cell_order(mesh.nCells(), -1);
    // Fill with consecutive numbers
    std::iota(cell_order.begin(), cell_order.end(), 0);

    // sort the cells in microdomain order
    std::stable_sort(cell_order.begin(), cell_order.end(), [&cellDist](int cell1, int cell2){
        return cellDist[cell1] < cellDist[cell2];
    });

    // Number of subdomains per MPI node
    int n_subdomains = n_domains / n_nodes;
    Foam::Info << "Number of microdomains per MPI node = " << n_subdomains << Foam::endl;
    if (n_domains % n_nodes != 0) {
        Foam::Info << "Error: Number of domains " << n_domains << " must be divisible by number of nodes " << n_nodes << Foam::endl;
        std::exit(1);
    }

    std::vector<int> node_order;
    for (int celli = 0; celli < mesh.nCells(); celli++) {
        int domain = int(cellDist[celli]);
        node       [celli] = domain / n_subdomains;
        microdomain[celli] = domain % n_subdomains;
        node_order.push_back(cellDist[cell_order[celli]] / n_subdomains);
    }

    node.write();
    microdomain.write();

    save_list("constant/cell_order.txt", cell_order);
    save_list("constant/node.txt",       node_order);
}
