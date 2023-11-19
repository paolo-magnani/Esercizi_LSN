#include "es10.h"
#define dim 1000
#define M 100000
#define N_migr 100
int main(int argc, char* argv[]){

    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int itag0 = 1, itag1 = 2, itag2 = 3, itag3 = 4;
    

    // i processi non comunicano
    cout << endl << size << " processi non comunicanti. " << endl;

    Pop_square pop(dim, "US_capitals.dat", rank); // ogni processo ha il generatore random inizializzato diversamente...
    //pop.print_city();
    pop.evolveL1(M);

    vector<double> best = pop.getchromo(0);
    cout << endl << "Processo rank: " << rank << " completato. Migliore distanza: " << best.back();


    string file1 = "indipendent/distance_rank_" + to_string(rank) + ".dat"; // ...e stampa su file differenti
    string file2 = "indipendent/path_rank_" + to_string(rank) + ".dat";
    pop.getresults(file1, file2);
    cout << endl << endl;

    // i processi comunicano
    cout << "Processi comunicanti ogni " << N_migr << " generazioni" << endl;

    Pop_square pop2(dim, "US_capitals.dat", rank);

    for(unsigned int i=0; i<N_migr; i++){

        MPI_Status stat1, stat2, stat3, stat4;
        MPI_Request req, req2;

        unsigned int n = M/N_migr;
        pop2.evolveL1(n);

        vector<double> v_send = pop2.getchromo(0);
        vector<double> v_receive(v_send.size());

        if(i%2==0){ // alternativamente, faccio comunicare tra loro due processi

            // scambio tra 1 & 0
            if(rank==1){
                // invio data
                MPI_Isend(v_send.data(), v_send.size(), MPI_DOUBLE_PRECISION, 0, itag0, MPI_COMM_WORLD, &req);
                //cout << " Sono " << rank <<  " e invio " << v_send.back() << endl;

                // aspetto di riceverla
                MPI_Recv(v_receive.data(), v_receive.size(), MPI_DOUBLE_PRECISION, 0, itag1, MPI_COMM_WORLD, &stat2);

            }else if(rank==0){

                MPI_Send(v_send.data(), v_send.size(), MPI_DOUBLE_PRECISION, 1, itag1, MPI_COMM_WORLD);
                //cout << " Sono " << rank <<  " e invio " << v_send.back() << endl;

                MPI_Recv(v_receive.data(), v_receive.size(), MPI_DOUBLE_PRECISION, 1, itag0, MPI_COMM_WORLD, &stat1);

            }

            // scambio tra 2 & 3
            if(rank==3){
                
                MPI_Isend(v_send.data(), v_send.size(), MPI_DOUBLE_PRECISION, 2, itag2, MPI_COMM_WORLD, &req2);
                //cout << " Sono " << rank <<  " e invio " << v_send.back() << endl;

                MPI_Recv(v_receive.data(), v_receive.size(), MPI_DOUBLE_PRECISION, 2, itag3, MPI_COMM_WORLD, &stat4);


            }else if(rank==2){

                MPI_Send(v_send.data(), v_send.size(), MPI_DOUBLE_PRECISION, 3, itag3, MPI_COMM_WORLD);
                //cout << " Sono " << rank <<  " e invio " << v_send.back() << endl;

                MPI_Recv(v_receive.data(), v_receive.size(), MPI_DOUBLE_PRECISION, 3, itag2, MPI_COMM_WORLD, &stat3);

            }

        }
        else{

            // scambio tra 1 & 3
            if(rank==1){

                MPI_Isend(v_send.data(), v_send.size(), MPI_DOUBLE_PRECISION, 3, itag0, MPI_COMM_WORLD, &req);
                //cout << " Sono " << rank <<  " e invio " << v_send.back() << endl;

                MPI_Recv(v_receive.data(), v_receive.size(), MPI_DOUBLE_PRECISION, 3, itag1, MPI_COMM_WORLD, &stat2);

            }else if(rank==3){

                MPI_Send(v_send.data(), v_send.size(), MPI_DOUBLE_PRECISION, 1, itag1, MPI_COMM_WORLD);
                //cout << " Sono " << rank <<  " e invio " << v_send.back() << endl;

                MPI_Recv(v_receive.data(), v_receive.size(), MPI_DOUBLE_PRECISION, 1, itag0, MPI_COMM_WORLD, &stat1);

            }

            // scambio tra 2 & 0
            if(rank==0){

                MPI_Isend(v_send.data(), v_send.size(), MPI_DOUBLE_PRECISION, 2, itag2, MPI_COMM_WORLD, &req2);
                //cout << " Sono " << rank <<  " e invio " << v_send.back() << endl;

                MPI_Recv(v_receive.data(), v_receive.size(), MPI_DOUBLE_PRECISION, 2, itag3, MPI_COMM_WORLD, &stat4);


            }else if(rank==2){

                MPI_Send(v_send.data(), v_send.size(), MPI_DOUBLE_PRECISION, 0, itag3, MPI_COMM_WORLD);
                //cout << " Sono " << rank <<  " e invio " << v_send.back() << endl;

                MPI_Recv(v_receive.data(), v_receive.size(), MPI_DOUBLE_PRECISION, 0, itag2, MPI_COMM_WORLD, &stat3);

            }

        }

        pop2.replacechromo(v_receive, dim-1);
        pop2.gensort();

    }

    cout << endl << "Processo rank: " << rank << " completato. Migliore distanza: " << best.back();
    file1 = "communicating/distance_rank_" + to_string(rank) + ".dat"; 
    file2 = "communicating/path_rank_" + to_string(rank) + ".dat";
    pop2.getresults(file1, file2);


    MPI_Finalize();
    return 0;
}