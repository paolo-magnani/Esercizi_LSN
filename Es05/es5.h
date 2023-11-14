# include "../libremia.h"

struct input{ //struct di input
	bool Gauss; // = 1 se si estraggono numeri da una distribuzione gaussiana, = 0 se si estraggono da una distribuzione uniforme 
	bool psi; // = 1 per lo stato 2p, = 0 per lo stato 1s
	double step; // step per avere l'accettanza a ~0.5
    double xi, yi, zi; // posizione iniziale
	string mean_data; // nome file di output
	string psi_data;
};

void ReadInput(input& myinput){ // funzione che legge i file di input

    ifstream in;
    in.open("input.in");
	if(!in.is_open()){
		cerr << "PROBLEM: unable to read input file" << endl;
	}else{
	in >> myinput.Gauss;
	in >> myinput.psi;
	in >> myinput.step;
	in >> myinput.xi >> myinput.yi >> myinput.zi;
	in >> myinput.mean_data;
	in >> myinput.psi_data;
	}
    in.close();
}

double radius(const double& x, const double& y, const double& z){ // restituisce il valore del raggio date le tre coordinate cartesiane
	return sqrt((x*x)+(y*y)+(z*z));
}

double psi100(const double& x, const double& y, const double& z){ // restituisce la funzione di probabilità dell'orbitale 1 0 0
	double psi = exp(-sqrt((x*x)+(y*y)+(z*z)))/(sqrt(M_PI));
	return psi;
}


double psi210(const double& x, const double& y, const double& z){ // restituisce la funzione di probabilità dell'orbitale 2 1 0
	double r = sqrt((x*x)+(y*y)+(z*z));
	if(r==0) return 0;
	
	double cos = z/r;
	double norm = sqrt(2./M_PI)/8.;
	double esp = exp(-r/2.);
	double psi = norm*r*esp*cos;
	return psi;
}


double accept(const double ratio){ // restituisce la probabilità di accettazione di una mossa
	if(ratio<1.) return ratio;
	else return 1.;
}