#include "es10.h"

int main(){

    MPI_Init(NULL,NULL);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat1, stat2;
    MPI_Request req;
    
    Pop_square pop(1000, "US_capitals.dat");
    pop.print_city();
    pop.evolveL1(1000000);
    pop.getresults();

    return 0;
}