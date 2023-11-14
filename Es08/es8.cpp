#include "es8.h"
#define M 100000
#define N 100

int main(int argc, char** argv){

	double step=1;

	if(argc<2)	cout << " <nstep> not specified, using default step instead" << endl;
	else step = atof(argv[1]);
	
	Random rand;
	iniz(rand);
	ofstream out_beta, out_par, out_histo;
	out_beta.open("Beta.dat");
	out_par.open("Parameters.dat");
	
	double ratio, weight_old, weight_new;
	double mu=rand.Rannyu(), sigma=rand.Rannyu();
	
	vector<double> E(M, 0.), mean(N, 0.), err(N, 0.);
	
	//Simulated Annealing
	for(double beta=1; beta<=1000; beta+=3){ //slowly freezing
		
		unsigned int Tstep=1000/sqrt(beta);
	
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

	stampamedia("Energy.dat",M, mean,err);

	out_beta.close();
	out_par.close();
	
	return 0;
}
