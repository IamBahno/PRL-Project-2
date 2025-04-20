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

// flatten vector of neighboars into an array of int
void flatten_vector(int *array, vector<Neighbour> *vector)
{
    int index = 0;
    for (int i = 0; i < vector->size(); i++)
    {
        array[index] = vector->at(i).forward;
        array[index + 1] = vector->at(i).reverse;
        
        index += 2;
    }
}
void send_neighbours(vector<vector<Neighbour>> *neighbours, int size){
    int neighbours_size = neighbours->size();
    // Broadcast the total number of vectors to all processes
    MPI_Bcast(&neighbours_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for(int vector = 0; vector < neighbours_size; vector++){
        int vector_size = neighbours->at(vector).size() *2;//how many ints are there in thsi vector;
        int *int_array = new int[vector_size];
        flatten_vector(int_array, &neighbours->at(vector));
        // Broadcast the size of the vector (number of elements)
        MPI_Bcast(&vector_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // Broadcast the flattened vector
        MPI_Bcast(int_array, vector_size, MPI_INT, 0, MPI_COMM_WORLD);
        
        delete[] int_array;
    }
}

void receive_neighbours(vector<vector<Neighbour>> *neighbours,int rank) {
    // Receive the total number of vectors (neighbours_size)
    int neighbours_size;
    MPI_Bcast(&neighbours_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Now, receive each vector and reconstruct it
    for (int vector = 0; vector < neighbours_size; vector++) {
        int vector_size;
        
        // Receive the size of the current vector
        MPI_Bcast(&vector_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Allocate memory for the flattened array
        int *int_array = new int[vector_size];

        // Receive the flattened vector
        MPI_Bcast(int_array, vector_size, MPI_INT, 0, MPI_COMM_WORLD);

        // Reconstruct the vector (convert int_array back into a vector of Neighbour)
        std::vector<Neighbour> current_vector;
        for (int i = 0; i < vector_size; i += 2) {
            // Create a Neighbour from pairs of integers
            Neighbour neighbour = { int_array[i], int_array[i + 1] };
            current_vector.push_back(neighbour);
        }

        // Add the reconstructed vector to the neighbours list
        neighbours->push_back(current_vector);

        // Clean up allocated memory for the flattened array
        delete[] int_array;
    }
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

    int edge_id;
    vector<vector<Neighbour>> new_neighbours;

    if (rank == 0) {
        create_adjecency_list(tree_string,&edges,&neighbours);

        send_neighbours(&neighbours,size);
        new_neighbours = neighbours;
    } 
    else{
        receive_neighbours(&new_neighbours,rank);
    }


    // Finalize
    MPI_Finalize();
    return 0;
}

// TODO find out how many procceses to be spawned