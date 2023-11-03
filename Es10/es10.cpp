#include "genalg.h"

int main(){

    MPI_Init(NULL,NULL);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat1, stat2;
    MPI_Request req;
    
    ifstream input;
    
    input.open("US_capitals.dat");
    cout<<"Carico le città dal file ..."<<endl;
    
    american_city cities;
    
    int n=0;
    while(!input.eof()){
        string state, name;
        double x,y;
        input>>x>>y;
        cities.x.push_back(x);
        cities.y.push_back(y);
        n++;
        //cout << endl << n << endl;
    }

    cout << endl << "Operazione completata. " << endl;

    Pop pop(200, cities);
    pop.print_city();

    input.close();
}