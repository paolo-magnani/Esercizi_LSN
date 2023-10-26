#include "../libremia.h"
#define M 1000000
#define N 100
//#define step 3

struct input{
	bool Gauss;
	bool psi;
	double step;
	string mean_data;
	string psi_data;
};

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



int main(){
	
	input myinput;
	ifstream in;
	in.open("input.in");
	if(!in.is_open()){
		cerr << "PROBLEM: unable to read input file" << endl;
	}else{
	in >> myinput.Gauss;
	in >> myinput.psi;
	in >> myinput.step;
	in >> myinput.mean_data;
	in >> myinput.psi_data;
	}
	
	const unsigned int L=M/N;
	Random rand;
	iniz(rand);
	ofstream outpsi, outmean;
	outpsi.open(myinput.psi_data);
	outmean.open(myinput.mean_data);
	
	double psi_old, psi_new, ratio, alpha;
	int acc=0; 
	//bool check;
	double xi = 0., yi = 0., zi = 0.;
	vector<double> x(M,0.), y(M,0.), z(M,0.), r(M,0.), media(N, 0.), sigma(N, 0.);
	
	
	for(unsigned int i=0; i<M; i++){
		
		if(myinput.psi==true){
			psi_old = pow(psi210(xi,yi,zi),2.);
		}else{
			psi_old = pow(psi100(xi,yi,zi),2.);
		}
		
		x[i]=xi;
		y[i]=yi;
		z[i]=zi;
		r[i]=radius(xi,yi,zi);
		
		if(myinput.Gauss==false){
			xi = xi + myinput.step*rand.Rannyu(-1.,1.);
			yi = yi + myinput.step*rand.Rannyu(-1.,1.);
			zi = zi + myinput.step*rand.Rannyu(-1.,1.);
		}else{
			xi = rand.Gauss(xi,myinput.step);
			yi = rand.Gauss(yi,myinput.step);
			zi = rand.Gauss(zi,myinput.step);
		}
		
		if(myinput.psi==true){
			psi_new = pow(psi210(xi,yi,zi),2.);
		}else{
			psi_new = pow(psi100(xi,yi,zi),2.);
		}
		
		
		if(psi_old == 0) ratio=1;
		else ratio = psi_new/psi_old;
		
		alpha = accept(ratio);
		if(rand.Rannyu()<=alpha){
			acc++;
			outpsi << xi << "," << yi << "," << zi << endl;
			}
		else{
			xi=x[i];
			yi=y[i];
			zi=z[i];
		}
		
	}
	mediablocchi(M, N, r, media, sigma);
	
	for(unsigned int i=0; i<N; i++) outmean << i*L << "," << media[i] << "," << sigma[i] << endl;
	cout << "Mean value: " << media[N-1] << endl;
	cout << "Acceptance rate: " << double(acc)/double(M) << endl;
	outpsi.close();
	outmean.close();
	in.close();
	
	return 0;
}
