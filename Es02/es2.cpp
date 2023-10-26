#include "../libremia.h"
#define M 10000000
#define N 100
#define ITER 100000 
#define NBLOC 100


void AcceptReject(vector<double> &, unsigned int);

void RandomWalk(bool, unsigned int, vector<double> &, Random  &); 

double error(const double &, const double &, unsigned int);

int main(){
	
	const unsigned int L=M/N;
	Random rand;
	
	vector<double> Int(M), sum_prog(N,0), err_prog(N,0);
	
	//Esercizio 1, punto 1: uniform sampling
	iniz(rand);
	
	ofstream out;
	out.open("unif.dat");
	
	generate(Int.begin(),Int.end(),[&](){return (M_PI/2.)*cos(M_PI*rand.Rannyu()/2.);});
	
	mediablocchi(M, N, Int, sum_prog, err_prog);
	
	for(unsigned int i=0; i<N; i++) out << (i+1.)*L << "," << sum_prog[i] << "," << err_prog[i] << endl;
	
	out.close(); 
	
	//Esercizio 1, punto 2: importance sampling
	
	ofstream out2;
	out2.open("impsamp.dat");
	
	vector<double> r(M, 0);
	
	AcceptReject(r, M);
	
	//for(unsigned int i=0; i<M; i++) Int[i] = (5.*M_PI/12.)*cos(M_PI*r[i]/2.)/(1.-(r[i]*r[i]/2.));
	
	for(unsigned int i=0; i<M; i++) Int[i] = (M_PI/4.)*cos(M_PI*r[i]/2.)/(1.-r[i]);
	
	mediablocchi(M, N, Int, sum_prog, err_prog);
	
	for(unsigned int i=0; i<N; i++) out2 << (i+1.)*L << "," << sum_prog[i] << "," << err_prog[i] << endl;
	
	out2.close();
	
	
	//Esercizio 2
	const int BLOC = ITER/NBLOC;
	vector<double> ave(N,0.), tot(N,0.), tot2(N,0.), err(N,0.);
	vector<double> dist(N,0.);
	
	//punto 1
	
	ofstream out3;
	out3.open("RWdiscr.dat");
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
		
		out3 << (i+1.) << "," << tot[i] << "," << err[i] << endl;	
	}
	
	out3.close();
	
	
	//Esercizio 2, punto 2
	ofstream out4;
	out4.open("RWcont.dat");
	
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
		
		out4 << (i+1.) << "," << tot[i] << "," << err[i] << endl;	
	}
	
	out4.close();
	
	 
	
	return 0;
	
}










void AcceptReject(vector<double> & r, unsigned int K){

	Random rand;
	iniz(rand);
		
	for(unsigned int i=0; i<K; i++){
		
		double x, y, f;
		
		do{
			x = rand.Rannyu();
			f = 2*(1-x);
			y = rand.Rannyu(0.,2.);
		
			r[i] = x;
			
		}while(y>f);
		
	}
}



void RandomWalk(bool cont, unsigned int K, vector<double>& dist, Random & rand){
	
	double x=0, y=0, z=0;
	
	if(cont==false){
		
		double sign;
		
		for(unsigned int i=0; i<K; i++){
			
			double r = rand.Rannyu(); 
			double r2 = rand.Rannyu();
			
			if(r2<0.5) sign = +1.;
			else sign = -1.;
			
			if(r<(1./3.)) x += sign;
			else if(r>(1./3.) and r<(2./3.)) y += sign;
			else z += sign;
			
			dist[i] = pow(x,2.)+pow(y,2.)+pow(z,2.);			
		}	
			
	
	}
	else{
		
		
		for(unsigned int i=0; i<K; i++){
		
			double cosine = rand.Rannyu(-1.,1.);
			double phi = rand.Rannyu(0., 2*M_PI);
			
			x += sqrt(1.-pow(cosine, 2.))*cos(phi);
			y += sqrt(1.-pow(cosine, 2.))*sin(phi);
			z += cosine;
			
			dist[i] = pow(x,2.)+pow(y,2.)+pow(z,2.);			
		}
		
		
	}
}


double error(const double & AV, const double & AV2, unsigned int n){
	if(n==0) return 0;
	else return sqrt((AV2 - pow(AV,2.))/double(n));
}

/*
void AcceptReject(vector<double> & r, unsigned int K){

	Random rand;
	iniz(rand);
		
	for(unsigned int i=0; i<K; i++){
		
		double x, y, f;
		
		do{
			x = rand.Rannyu();
			f = (6./5.)*(1.-(x*x/2.));
			y = rand.Rannyu(0.,(6./5.));
		
			r[i] = x;
			
		}while(y>f);
		
	}
}*/
