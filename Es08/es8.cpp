#include "../libremia.h"
#define M 100000
#define N 100


double psiplus(const double & x, const double & mu, const double & sigma){
	return exp(-0.5*pow((x+mu)/sigma, 2.));
}

double psiminus(const double & x, const double & mu, const double & sigma){
	return exp(-0.5*pow((x-mu)/sigma, 2.));
}

double psi(const double & x, const double & mu, const double & sigma){
	return psiminus(x,mu,sigma)+psiplus(x,mu,sigma);;
}


double Ham(const double & x, const double & mu, const double & sigma){
	double V = (((x*x)-2.5)*x*x);
	double K = -0.5*(mu*mu + (x*x)- (sigma*sigma)+ (2*mu*x*(psiplus(x, mu, sigma)-psiminus(x, mu, sigma))/psi(x, mu, sigma)))*pow(1./sigma,4.);
	return V+K;
}

double Boltz(const double & E, const double & beta){
	return exp(-beta*E);
	
} 

double accept(const double ratio){
	if(ratio<1.) return ratio;
	else return 1.;
}


double Eval_H(const unsigned int& B, const double& step, const double& mu, const double& sigma, Random& rand){

	double E=0;
	double xi=0., psi_old, psi_new, ratio;
	
	for(unsigned int i=0; i<B; i++){
			
			psi_old = pow(psi(xi,mu,sigma),2.);
			
			double x_old = xi; 
			
			E += Ham(xi, mu, sigma);
			
			xi = xi + step*rand.Rannyu(-1.,1.);
			
			psi_new = pow(psi(xi,mu,sigma),2.);
			
			if(psi_old == 0) ratio=1.;
			else ratio = psi_new/psi_old;
			
			double alpha = accept(ratio);
			if(rand.Rannyu()>alpha) xi=x_old;
	
	}

	return E/double(B);

} 


void Stima_Finale(const unsigned int& B, double& step, const double& mu, const double& sigma, vector<double>& E, Random& rand){

	double xi=0., psi_old, psi_new, ratio;
	double acc = 0, acc_old;

	do{
		acc_old = acc;
		acc = 0, xi=0;	
		for(unsigned int i=0; i<100; i++){
			
			psi_old = pow(psi(xi,mu,sigma),2.);
			
			double x_old = xi;			
			xi = xi + step*rand.Rannyu(-1.,1.);
			
			psi_new = pow(psi(xi,mu,sigma),2.);
			
			if(psi_old == 0) ratio=1.;
			else ratio = psi_new/psi_old;
			
			double alpha = accept(ratio);
			if(rand.Rannyu()<=alpha) acc++;
			else xi=x_old;
	
		}
		acc /= double(100);
		if(acc<acc_old) step -= 0.1;
		else step += 0.1;
		cout << "acceptance: " << acc << endl;
		cout << "step: " << step << endl;
	}while(acc<0.4 or acc>0.6);
	
	acc = 0;
	for(unsigned int i=0; i<B; i++){
			
			psi_old = pow(psi(xi,mu,sigma),2.);
			
			double x_old = xi; 
			
			E[i]= Ham(xi, mu, sigma);
			
			xi = xi + step*rand.Rannyu(-1.,1.);
			
			psi_new = pow(psi(xi,mu,sigma),2.);
			
			if(psi_old == 0) ratio=1.;
			else ratio = psi_new/psi_old;
			
			double alpha = accept(ratio);
			if(rand.Rannyu()<=alpha) acc++;
			else xi=x_old;
	
	}
	acc /= double(B);
}



int main(int argc, char** argv){
	
	const unsigned int L=M/N;
	double step=1;

	if(argc<2)	cout << " <nstep> not specified, using default step instead" << endl;
	else step = atof(argv[1]);
	
	Random rand;
	iniz(rand);
	ofstream out_beta, out_par, out_H, out_histo;
	out_beta.open("Beta.dat");
	out_par.open("Parameters.dat");
	out_H.open("Energy.dat");
	
	double psi_old, psi_new, ratio, weight_old, weight_new;
	double mu=rand.Rannyu(), sigma=rand.Rannyu();
	double xi=0.;
	
	vector<double> E(M, 0.), mean(N, 0.), err(N, 0.);
	
	//Simulated Annealing
	for(double beta=1; beta<=1000; beta+=3){ //slowly freezing
		
		unsigned int Tstep=1000/beta;
	
		//MU, SIGMA
		for(unsigned int i=0; i<Tstep; i++){
			
			cout << " beta = " << beta << "; Avanzamento:" << (double(i+1)/Tstep)*100 << "%" << endl;

			//HAMILTONIAN VALUE
			double H_old = Eval_H(M, step, mu, sigma, rand);
			weight_old = Boltz(H_old, beta);
			
			double mu_old = mu;
			double sigma_old = sigma;
			
			mu = fabs(mu + (1/beta)*rand.Rannyu(-1,1));
			sigma = fabs(sigma + (1/beta)*rand.Rannyu(-1,1));
			
			double H_new = Eval_H(M, step, mu, sigma, rand);;
			
			weight_new = Boltz(H_new, beta);
			
			if(weight_old == 0) ratio=1.;
				else ratio = weight_new/weight_old;
				
				double alpha = accept(ratio);
				if(rand.Rannyu()>alpha){
					mu=mu_old;
					sigma = sigma_old;
				}	
		}
		double H =  Eval_H(M, step, mu, sigma, rand);

		out_par << beta << "," << mu << "," << sigma << endl;
		out_beta << beta << "," << H << endl;

		cout << "Mu: " << mu << endl;
		cout << "Sigma: " << sigma << endl;
		cout << "Estimated < H > value: " << H << endl;
	}
	Stima_Finale(M, step, mu, sigma, E, rand);
	mediablocchi(M, N, E, mean, err);

	for(unsigned int i=0; i<N; i++) out_H << i+1 << "," << mean[i] << "," << err[i] << endl;

	out_beta.close();
	out_par.close();
	out_H.close();
	
	return 0;
}
