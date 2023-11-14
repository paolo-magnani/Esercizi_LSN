#include "es5.h"
#define M 1000000 // numero di step per il Metropolis
#define N 100 // numero di blocchi


int main(){

	input myinput;

	ReadInput(myinput); // prendo variabili di input dal file input.in

	Random rand;
	iniz(rand); // inizializzo il generatore di numeri casuali
	ofstream outpsi;
	outpsi.open(myinput.psi_data);
	
	double psi_old, psi_new, ratio, alpha;
	int acc=0; // inizializzo l'accettanza a 0
	double xi=myinput.xi, yi=myinput.yi, zi=myinput.zi;
	vector<double> x(M,0.), y(M,0.), z(M,0.), r(M,0.), media(N, 0.), sigma(N, 0.);
	
	
	for(unsigned int i=0; i<M; i++){ // implementazione dell'algoritmo di Metropolis
		
		if(myinput.psi==true){
			psi_old = pow(psi210(xi,yi,zi),2.); // densità di probabilità del punto di partenza
		}else{
			psi_old = pow(psi100(xi,yi,zi),2.);
		}
		
		x[i]=xi;
		y[i]=yi;
		z[i]=zi;
		r[i]=radius(xi,yi,zi); // registra il valore della posizione
		
		if(myinput.Gauss==false){ // tento una mossa in una generica direzione e secondo una distribuzione definita (uniforme o gaussiana)
			xi = xi + myinput.step*rand.Rannyu(-1.,1.);
			yi = yi + myinput.step*rand.Rannyu(-1.,1.);
			zi = zi + myinput.step*rand.Rannyu(-1.,1.);
		}else{
			xi = rand.Gauss(xi,myinput.step);
			yi = rand.Gauss(yi,myinput.step);
			zi = rand.Gauss(zi,myinput.step);
		}
		
		if(myinput.psi==true){ // dalla nuova posizione calcolo la distribuzione di probabilità
			psi_new = pow(psi210(xi,yi,zi),2.);
		}else{
			psi_new = pow(psi100(xi,yi,zi),2.);
		}
		
		
		if(psi_old == 0) ratio=1; // per evitare errori, la mossa è sempre accettata quando la posizione vecchia dà probabilità = 0
		else ratio = psi_new/psi_old;
		
		alpha = accept(ratio);
		if(rand.Rannyu()<=alpha){ // accetta la mossa secondo una probabilità che dipende la rapporto tra le densità di probabilità
			acc++;
			outpsi << xi << "," << yi << "," << zi << endl;
			}
		else{ // altrimenti rimane nella stessa posizione
			xi=x[i];
			yi=y[i];
			zi=z[i];
		}
		
	}
	mediablocchi(M, N, r, media, sigma); // calcolo il valor medio del raggio

	stampamedia(myinput.mean_data, M, media, sigma); // stampo dati su file

	cout << "Mean value: " << media[N-1] << endl;
	cout << "Acceptance rate: " << double(acc)/double(M) << endl; // stampo a video il valore dell'accettanza
	outpsi.close();

	
	return 0;
}
