#include "es10.h"
#define M 10000
#define N_migr 100
int main(int argc, char* argv[]){

    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat1, stat2;
    MPI_Request req;
    

    // i processi non comunicano
    Pop_square pop(1000, "US_capitals.dat", rank); // ogni processo ha il generatore random inizializzato diversamente...
    //pop.print_city();
    pop.evolveL1(10000);
    string file1 = "indipendent/distance_rank_" + to_string(rank) + ".dat"; // ...e stampa su file differenti
    string file2 = "indipendent/path_rank_" + to_string(rank) + ".dat";
    pop.getresults(file1, file2);

    // i processi comunicano
    // Pop_square pop2(1000, "US_capitals.dat", rank);

    // for(unsigned int i=0; i<N_migr; i++){
    //     unsigned int n = M/N_migr;
    //     pop2.evolveL1(n);

    //     if

    // }


    MPI_Finalize();
    return 0;
}