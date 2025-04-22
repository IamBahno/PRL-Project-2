#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>

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

// flatten vector of neighbours into an array of int
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

// broadcast the adjecency list to all processes
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

// Receive and reconstruct the adjecency list
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

// Printing for debugging
void print_adjecency_list(vector<vector<Neighbour>> new_neighbours){
    for (size_t i = 0; i < new_neighbours.size(); ++i) {
        std::cerr << "["; 
        for (size_t j = 0; j < new_neighbours[i].size(); ++j) {
            const auto& neighbour = new_neighbours[i][j];
            
            std::cerr << neighbour.forward << "," << neighbour.reverse;

            if (j != new_neighbours[i].size() - 1) {
                std::cerr << ", ";
            }
        }
        std::cerr << "] ";
    }
    std::cerr << std::endl;

}


// Returns the next edge to be visited  in euler tour
int create_euler_tour(int edge, vector<vector<Neighbour>>neighbours){
    // First find the reverse edge to the currente edge
    int reverse;
    for(int i =0; i < neighbours.size();i++){
        for(int j = 0; j < neighbours.at(i).size();j++){
            if(edge == neighbours.at(i).at(j).forward){
                reverse = neighbours.at(i).at(j).reverse;
            }
        }
    }

    // Now find where is the reverse edge as forward edge
    for(int i =0; i < neighbours.size();i++){
        for(int j = 0; j < neighbours.at(i).size();j++){
            if(reverse == neighbours.at(i).at(j).forward){
                // If its not last last edge pair at given list, return the next forward edge after the found pair
                if(j < neighbours.at(i).size() - 1){
                    return neighbours.at(i).at(j+1).forward;
                }
                // Else return the the first forwas edge at given edge pair list
                else{
                    return neighbours.at(i).at(0).forward;
                }
            }
        }
    }
    return -1;
}

// Rank 0 collect the euler root into *euler_toor
void euler_to_rank0(vector<int> *euler_tour,int next_tour,int rank,int world_size){
    if(rank == 0){
        euler_tour->at(0) = next_tour; // sets its own value
        int next_edge;
        for(int i = 1; i < world_size;i++)
        {
            MPI_Recv(&next_edge, 1, MPI_INT, i, 0, MPI_COMM_WORLD, NULL);
            euler_tour->at(i) = next_edge;
        }
    }
    else{
        MPI_Send(&next_tour, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

// The last edge that is going to the root (the second one) is looped to the same node it coming from
// that way creating a root
// Return the id of the self loop edge
int introduce_root(vector<int> *euler_tour,vector<Edge> *edges,char root_char){
    bool first = true;
    for(Edge& edge : *edges){
        if(edge.to == root_char){
            if(first == true){first = false;continue;}
            edge.to = edge.from;
            euler_tour->at(edge.id) = edge.id;
            return edge.id;
        }
    }
    return -1;
}

// ranks is the vector which will be filled with the distance to end (ranks)
// parameter rank is id of process/edge
// euler tour is euler tour with root altered
// Using the path doubling alghorithm
void create_rank_vector(vector<int> *ranks,vector<int> *euler_tour,int rank,int size, int self_loop_edge_id){
    
    //crete copy of euler rout
    //we dont want cange the original
    vector<int> succesor_list;
    if (rank == 0) {
        succesor_list = *euler_tour;
    }else{
        succesor_list.resize(size);
    }
    // share it with all processes
    MPI_Bcast(succesor_list.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
    // rank 0 creates a ranks vector, filled with ones, only the self loop edge will be 0
    if (rank == 0) {
        ranks->assign(size, 1);
        ranks->at(self_loop_edge_id) = 0;
    }else{
        ranks->resize(size);
    }
    // share with all
    MPI_Bcast(ranks->data(), size, MPI_INT, 0, MPI_COMM_WORLD);

    // each edge in each iteration find its rank and succesor
    // all processes then then gather the new values together, so each process have new updated vectors 
    for (int k = 1; k < ceil(log2(size)) + 1;k++){
        int new_rank = ranks->at(rank) + ranks->at(succesor_list.at(rank));
        int new_succ = succesor_list.at(succesor_list.at(rank));
        MPI_Allgather(&new_rank, 1, MPI_INT, ranks->data(), 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&new_succ, 1, MPI_INT, succesor_list.data(), 1, MPI_INT, MPI_COMM_WORLD);
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

    vector<vector<Neighbour>> new_neighbours;

    if (rank == 0) {
        // Rank 0 parses the input string and creates the adjacency list
        // I use a vector for the adjacency list since it's conviniet dynamic type
        create_adjecency_list(tree_string,&edges,&neighbours);

        //  and distributes neighbour information to all processes
        send_neighbours(&neighbours,size);
        new_neighbours = neighbours;
    } 
    else{
        // Other proccesses receives the neighbour list
        receive_neighbours(&new_neighbours,rank);
    }

    //Each process/edge finds its following edge 
    int edge_id = rank;
    int next_tour = create_euler_tour(edge_id,new_neighbours);

    // Rank 0 collects the euler tour to euler_tour
    vector<int> euler_tour(edges.size());
    euler_to_rank0(&euler_tour,next_tour,rank,size);

    int self_loop_edge_id;
    //Change the root so there is root
    if(rank == 0){
        self_loop_edge_id = introduce_root(&euler_tour,&edges,tree_string[0]);
    }
    // Send the euler tour values backto the proccesse, one proccess will recieve changed value
    MPI_Scatter(euler_tour.data(),1,MPI_INT,&next_tour,1,MPI_INT,0,MPI_COMM_WORLD);

    // Compute rank vector (rank is the distance of each edge to the end)
    vector<int> ranks;
    create_rank_vector(&ranks,&euler_tour,rank,size,self_loop_edge_id);

    
    // if(rank == 0){
    //     for(int i = 0; i < size;i++){
    //         std::cout << ranks.at(i)<< std::endl;
    //     }
    // }
    // if(rank == 0){
    //     for(int i = 0; i < euler_tour.size();i++){
    //         std::cout << "Edge: "<< i << "next: " << euler_tour.at(i)<< std::endl;
    //     }
    // }
    // std::cerr << "Rank " << rank << " received euler_tour value: " << next_tour << std::endl;


    // Finalize
    MPI_Finalize();
    return 0;
}