#include "../libremia.h"
#define M 100000
#define N 100



double psi(const double & x, const double & mu, const double & sigma){
	double func = exp(-0.5*pow((x-mu)/sigma, 2.))+exp(-0.5*pow((x+mu)/sigma, 2.));
	return func;
}

double Ham(const double & x, const double & mu, const double & sigma, const double & psi){

	//double E = (((x*x)-2.5)*x*x) - 0.5*(pow(1./sigma,2.)*(-1. + pow(x/sigma, 2.) + pow(mu/sigma, 2.) + (2*mu*x/(sigma*sigma))*(1. - (2.*exp(-0.5*pow((x-mu)/sigma, 2.))/psi))));
	double E = (x*x) - 0.5*(pow(1./sigma,2.)*(-1. + pow(x/sigma, 2.) + pow(mu/sigma, 2.) + (2*mu*x/(sigma*sigma))*(1. - (2.*exp(-0.5*pow((x-mu)/sigma, 2.))/psi)))); //potenziale test
	return E;
}

double Boltz(const double & E, const double & T){
	double b = 1./T;
	return exp(-b*E);
	
} 

double accept(const double ratio){
	if(ratio<1.) return ratio;
	else return 1.;
}


void Eval_H(const unsigned int B, const double& step, const double& mu, const double& sigma, vector<double>& E, Random& rand, int& acc){

	double xi=0., psi_old, psi_new, ratio;
	
	for(unsigned int i=0; i<B; i++){
			
			psi_old = pow(psi(xi,mu,sigma),2.);
			
			double x_old = xi; 
			
			E[i]= Ham(xi, mu, sigma, psi_old);
			
			
			xi = xi + step*rand.Rannyu(-1.,1.);
			
			psi_new = pow(psi(xi,mu,sigma),2.);
			
			if(psi_old == 0) ratio=1.;
			else ratio = psi_new/psi_old;
			
			double alpha = accept(ratio);
			if(rand.Rannyu()<=alpha) acc++;
			else xi=x_old;
	
	}

} 


int main(int argc, char** argv){

	double step=2.1;

	if(argc<2)	cout << " <nstep> not specified, using default step instead" << endl;
	else step = atof(argv[1]);
	
	const unsigned int L=M/N;
	Random rand;
	iniz(rand);
	ofstream out;
	out.open("test.dat");
	
	double psi_old, psi_new, ratio, weight_old, weight_new;
	int acc=0;
	double mu=rand.Rannyu(), sigma=rand.Rannyu();
	double xi=0.;
	
	vector<double> E(M, 0.), mean(N, 0.), err(N, 0.);
	
	//BETA
	cout << endl << "check" << endl;
	
	for(double T=10; T>0; T=T-0.1){ //slowly freezing
		cout << " T = " << T;
		cout << endl << " Simulation is at " << T << "% ";
		unsigned int Tstep=1000*T;
		acc = 0;
	
		//MU, SIGMA
		for(unsigned int i=0; i<Tstep; i++){
			
			cout << endl << " Tstep = " << i << endl;
			//HAMILTONIAN VALUE
			for(unsigned int i=0; i<M; i++){
			
				psi_old = pow(psi(xi,mu,sigma),2.);
				
				double x_old = xi; 
				
				E[i]= Ham(xi, mu, sigma, psi_old);
				
				
				xi = xi + step*rand.Rannyu(-1.,1.);
				
				psi_new = pow(psi(xi,mu,sigma),2.);
				
				if(psi_old == 0) ratio=1.;
				else ratio = psi_new/psi_old;
				
				double alpha = accept(ratio);
				if(rand.Rannyu()<=alpha) acc++;
				else xi=x_old;
		
			}
			mediablocchi(M, N, E, mean, err);
			
			double H_old = mean[N-1];
			
			weight_old = Boltz(H_old, T);
			
			double mu_old = mu;
			double sigma_old = sigma;
			
			mu += (T/10)*rand.Rannyu(-1,1);
			sigma += (T/10)*rand.Rannyu(-1,1);
	
	
			for(unsigned int i=0; i<M; i++){
			
				psi_old = pow(psi(xi,mu,sigma),2.);
				
				double x_old = xi; 
				
				E[i]= Ham(xi, mu, sigma, psi_old);
				
				
				xi = xi + step*rand.Rannyu(-1.,1.);
				
				psi_new = pow(psi(xi,mu,sigma),2.);
				
				if(psi_old == 0) ratio=1.;
				else ratio = psi_new/psi_old;
				
				double alpha = accept(ratio);
				if(rand.Rannyu()<=alpha) acc++;
				else xi=x_old;
		
			}
			mediablocchi(M, N, E, mean, err);
			
			double H_new = mean[N-1];
			
			weight_new = Boltz(H_new, T);
			
			if(weight_old == 0) ratio=1.;
				else ratio = weight_new/weight_old;
				
				double alpha = accept(ratio);
				if(rand.Rannyu()>alpha){
					mu=mu_old;
					sigma = sigma_old;
				}	
		}
		
		//for(unsigned int i=0; i<N; i++) out << (i*L)+L << "," << mean[i] << "," << err[i] << endl;
		cout << "Mu: " << mu << endl;
		cout << "Sigma: " << sigma << endl;
		cout << "Mean value: " << mean[N-1] << endl;
		cout << "Acceptance rate: " << double(acc)/double(M) << endl;
	}
	
	out.close();
	
	return 0;
}
