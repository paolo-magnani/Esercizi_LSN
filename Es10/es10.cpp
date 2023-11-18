#include "es10.h"

int main(int argc, char* argv[]){

    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat1, stat2;
    MPI_Request req;
    
    Pop_square pop(1000, "US_capitals.dat", rank);
    //pop.print_city();
    pop.evolveL1(10000);
    string file1 = "indipendent/distance_rank_" + to_string(rank) + ".dat";
    string file2 = "indipendent/path_rank_" + to_string(rank) + ".dat";
    pop.getresults(file1, file2);

    MPI_Finalize();
    return 0;
}