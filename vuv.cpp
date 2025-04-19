#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>

using namespace std;

struct Edge {
    int id;
    char from;
    char to;
};

// Neighbour of and tree node, is represented like edge to the neighbour and reverse edge
// Storing only ids of the edges
struct Neighbour{
    int forward;
    int reverse;
};

void create_adjecency_list(string tree_string,vector<Edge> *edges, vector<vector<Neighbour>> *neighbours){
    return;
}

int main(int argc, char** argv) {
    // Init MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector<Edge> edges;
    vector<vector<Neighbour>> neighbours;

    std::string tree_string = argv[1];

    if (rank == 0) {
        std::cout << "🌲 trunk" << tree_string << "\n";

        // TODO creates an adjecent list
        create_adjecency_list(tree_string,&edges,&neighbours);
        // and share with other proccesses
    }



    // Finalize
    MPI_Finalize();
    return 0;
}

// TODO find out how many procceses to be spawned