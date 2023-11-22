#include "../libremia.h"


double psiplus(const double & x, const double & mu, const double & sigma){ // parte della funzione di prova del ground state
	return exp(-0.5*pow((x+mu)/sigma, 2.));
}

double psiminus(const double & x, const double & mu, const double & sigma){ // parte della funzione di prova del ground state
	return exp(-0.5*pow((x-mu)/sigma, 2.));
}

double psi(const double & x, const double & mu, const double & sigma){ // funzione di prova per il ground state
	return psiminus(x,mu,sigma)+psiplus(x,mu,sigma);;
}


double Ham(const double & x, const double & mu, const double & sigma){ // calcola l'Hamiltoniana del sistema
	double V = (((x*x)-2.5)*x*x);
	double K = -0.5*(mu*mu + (x*x)- (sigma*sigma)+ (2*mu*x*(psiplus(x, mu, sigma)-psiminus(x, mu, sigma))/psi(x, mu, sigma)))*pow(1./sigma,4.);
	return V+K;
}

double Boltz(const double & E, const double & beta){ // calcola il peso di Boltzmann
	return exp(-beta*E);
	
} 

double accept(const double ratio){ // restituisce la probabilità di accettazione di una mossa del Metropolis
	if(ratio<1.) return ratio;
	else return 1.;
}


double Eval_H(const unsigned int& B, const double& step, const double& mu, const double& sigma, Random& rand){ // stima il valor medio di H 

	double E=0;
	double xi=0., psi_old, psi_new, ratio;
	
	for(unsigned int i=0; i<B; i++){ // Algoritmo di Metropolis
			
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


void Stima_Finale(const unsigned int& B, double& step, const double& mu, const double& sigma, vector<double>& E, vector<double>& pos, Random& rand){ // stima il valore di H nel ground state dopo aver trovato il valore di mu e sigma

	double xi=0., psi_old, psi_new, ratio;
	double acc = 0, acc_old;

	cout << endl << "Searching for a suitable step... " << endl;
	do{ // trova il corretto valore dello step per avere un'accettanza di ~0.5 
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
	}while(acc<0.45 or acc>0.55);
	
	cout << endl << " Evaluating < H >... " << endl;
	acc = 0; // una volta trovato il valore dello step,
	for(unsigned int i=0; i<B; i++){ //campiona il valor medio di H con l'algoritmo di Metropolis e registra le posizioni esplorate in 'pos'
			
			psi_old = pow(psi(xi,mu,sigma),2.);
			
			double x_old = xi; 
			pos.push_back(x_old);
			
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
