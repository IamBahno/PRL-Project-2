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
    int n = tree_string.size();
    neighbours->resize(n);

    int edge_id = 0;

    // in each iteratoin create edges with the nodes childrens, forward and reverse
    for (int i = 0; i < n; ++i) {
        int left = 2 * i + 1;
        int right = 2 * i + 2;

        if (left < n) {
            edges->push_back({edge_id, tree_string[i], tree_string[left]});
            neighbours->at(i).push_back({edge_id, edge_id + 1});

            edges->push_back({edge_id + 1, tree_string[left], tree_string[i]});
            neighbours->at(left).push_back({edge_id + 1, edge_id});

            edge_id += 2;
        }

        if (right < n) {
            edges->push_back({edge_id, tree_string[i], tree_string[right]});
            neighbours->at(i).push_back({edge_id, edge_id + 1});

            edges->push_back({edge_id + 1, tree_string[right], tree_string[i]});
            neighbours->at(right).push_back({edge_id + 1, edge_id});

            edge_id += 2;
        }
    }

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

        create_adjecency_list(tree_string,&edges,&neighbours);
    }


    // Finalize
    MPI_Finalize();
    return 0;
}

// TODO find out how many procceses to be spawned