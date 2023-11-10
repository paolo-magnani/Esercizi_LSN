# include "../libremia.h"

struct input{
	bool Gauss;
	bool psi;
	double step;
    double xi, yi, zi;
	string mean_data;
	string psi_data;
};

void ReadInput(input& myinput){

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

double radius(const double& x, const double& y, const double& z){
	return sqrt((x*x)+(y*y)+(z*z));
}

double psi100(const double& x, const double& y, const double& z){
	double psi = exp(-sqrt((x*x)+(y*y)+(z*z)))/(sqrt(M_PI));
	return psi;
}


double psi210(const double& x, const double& y, const double& z){
	double r = sqrt((x*x)+(y*y)+(z*z));
	if(r==0) return 0;
	
	double cos = z/r;
	double norm = sqrt(2./M_PI)/8.;
	double esp = exp(-r/2.);
	double psi = norm*r*esp*cos;
	return psi;
}


double accept(const double ratio){
	if(ratio<1.) return ratio;
	else return 1.;
}