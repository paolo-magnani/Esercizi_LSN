#include "../libremia.h"
#define M 1000000 // numero ripetizioni
#define N 100 // numero blocchi
#define nintervals 100 // numerp di intervalli per il caso discreto


int main(){
	
	// parametri
	const double S_0 = 100;
	const double T = 1;
	const double t = T/double(nintervals);
	const double K = 100;
	const double r = 0.1;
	const double sigma = 0.25;
	
	const unsigned int L=M/N;
	Random rand;
	iniz(rand); // inizializzo il generatore di numeri casuali
	
	
	ofstream outP;
	outP.open("risultati/3_direct_put.dat");
	
	ofstream outC;
	outC.open("risultati/3_direct_call.dat");
	
	ofstream outP2;
	outP2.open("risultati/3_discret_put.dat");
	
	ofstream outC2;
	outC2.open("risultati/3_discret_call.dat");
	
	vector<double> P(M), C(M), P2(M), C2(M); //preparo le variabili
	vector<double> call(N,0), err_call(N,0);
	vector<double> put(N,0), err_put(N,0);
	vector<double> call2(N,0), err_call2(N,0);
	vector<double> put2(N,0), err_put2(N,0);
	
	cout << " Calcolo la stima diretta e discretizzata di put e call..." << endl;
	for(unsigned int i=0; i<M; i++){

		//stima diretta 
		P[i]=0;
		C[i] = exp(-r*T)*(S_0*exp((r-(0.5*pow(sigma,2)))*T+(sigma*rand.Gauss(0, sqrt(T))))-K); //implemento l'algoritmo visto a lezione: valuto il prezzo e lo sottraggo al prezzo iniziale
			
		if(C[i]<0){ // se il prezzo al tempo T è inferiore al prezzo iniziale, scambia put con call
		
			P[i]=-C[i];
			C[i]=0;
		}

		//stima discretizzata
		double S = S_0;
		P2[i]=0;
		for(unsigned int j=1; j<=nintervals; j++) S = S*exp((r-(0.5*pow(sigma,2)))*t+(sigma*rand.Gauss(0,1)*sqrt(t))); //valuto S per ogni step discretizzato
		C2[i]= exp(-r*T)*(S-K);
		
		if(C2[i]<0){

			P2[i]=-C2[i];
			C2[i]=0;
		}	

	}
		
	mediablocchi(M, N, C, call, err_call); // calcolo le medie
	mediablocchi(M, N, P, put, err_put);
	mediablocchi(M, N, C2, call2, err_call2);
	mediablocchi(M, N, P2, put2, err_put2);
	
	
	for(unsigned int i=0; i<N; i++){ // stampo i risultati
		outC << (i+1)*L << "," << call[i] << "," << err_call[i] << endl; 
		outP << (i+1)*L << "," << put[i] << "," << err_put[i] << endl; 
		outC2 << (i+1)*L << "," << call2[i] << "," << err_call2[i] << endl; 
		outP2 << (i+1)*L << "," << put2[i] << "," << err_put2[i] << endl;
	}
	
	outP.close();
	outC.close();
	outP2.close();
	outC2.close();

	cout << " Completo." << endl;
	
 	return 0;
}
