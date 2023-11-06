#include "es2.h"
#define M 10000000
#define N 100
#define ITER 100000 
#define NBLOC 100


int main(){
	
	Random rand;
	
	vector<double> Int(M), sum_prog(N,0), err_prog(N,0);
	
	//Esercizio 1, punto 1: uniform sampling
	cout << " Risolvo il punto 1.1..." << endl;

	iniz(rand);
	
	generate(Int.begin(),Int.end(),[&](){return (M_PI/2.)*cos(M_PI*rand.Rannyu()/2.);});
	
	mediablocchi(M, N, Int, sum_prog, err_prog);

	stampamedia("unif.dat",M,sum_prog, err_prog);
	
	//Esercizio 1, punto 2: importance sampling	
	cout << " Risolvo il punto 1.2..." << endl;

	vector<double> r(M, 0);
	
	AcceptReject(r, M);
	
	for(unsigned int i=0; i<M; i++) Int[i] = (M_PI/4.)*cos(M_PI*r[i]/2.)/(1.-r[i]);
	
	mediablocchi(M, N, Int, sum_prog, err_prog);

	stampamedia("impsamp.dat", M, sum_prog, err_prog);
	
	
	//Esercizio 2
	const int BLOC = ITER/NBLOC;
	vector<double> ave(N,0.), tot(N,0.), tot2(N,0.), err(N,0.);
	vector<double> dist(N,0.);
	
	//punto 1
	cout << " Risolvo il punto 2.1..." << endl;
	
	ofstream out;
	out.open("RWdiscr.dat");
	iniz(rand);	
	
	
	for(unsigned int i=0; i<NBLOC; i++){ 
	
		for(unsigned int k=0; k<BLOC; k++){
	
			RandomWalk(false, N, dist, rand);
		
			for(unsigned int j=0; j<N; j++) ave[j] += dist[j];
			
		}
		
		for(unsigned int j=0; j<N; j++){
		
			ave[j] = sqrt(ave[j]/double(BLOC));
			tot[j] += ave[j];
			tot2[j] += ave[j]*ave[j];
		}
	
		fill(ave.begin(),ave.end(), 0);
		
		
	}
	
	for(unsigned int i=0; i<N; i++){
		
		tot[i] /= double(NBLOC);
		tot2[i] /= double(NBLOC);
		
		err[i] = error(tot[i], tot2[i], NBLOC-1);
		
		out << (i+1.) << "," << tot[i] << "," << err[i] << endl;	
	}
	
	out.close();
	
	
	//Esercizio 2, punto 2
	cout << " Risolvo il punto 2.2..." << endl;

	out.open("RWcont.dat");
	
	fill(tot.begin(), tot.end(), 0.);
	fill(tot2.begin(), tot2.end(), 0.);
	fill(err.begin(), err.end(), 0.);
	
	for(unsigned int i=0; i<NBLOC; i++){ 
	
		for(unsigned int k=0; k<BLOC; k++){
	
			RandomWalk(true, N, dist, rand);
		
			for(unsigned int j=0; j<N; j++) ave[j] += dist[j];
			
		}
		
		for(unsigned int j=0; j<N; j++){
		
			ave[j] = sqrt(ave[j]/double(BLOC));
			tot[j] += ave[j];
			tot2[j] += ave[j]*ave[j];
		}
	
		fill(ave.begin(),ave.end(), 0);
		
		
	}
	
	for(unsigned int i=0; i<N; i++){
		
		tot[i] /= double(NBLOC);
		tot2[i] /= double(NBLOC);
		
		err[i] = error(tot[i], tot2[i], NBLOC-1);
		
		out << (i+1.) << "," << tot[i] << "," << err[i] << endl;	
	}
	
	out.close();
	
	return 0;
	
}












