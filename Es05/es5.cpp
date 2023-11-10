#include "es5.h"
#define M 1000000
#define N 100


int main(){

	input myinput;

	ReadInput(myinput);

	Random rand;
	iniz(rand);
	ofstream outpsi;
	outpsi.open(myinput.psi_data);
	
	double psi_old, psi_new, ratio, alpha;
	int acc=0; 
	double xi=myinput.xi, yi=myinput.yi, zi=myinput.zi;
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

	stampamedia(myinput.mean_data, M, media, sigma);

	cout << "Mean value: " << media[N-1] << endl;
	cout << "Acceptance rate: " << double(acc)/double(M) << endl;
	outpsi.close();

	
	return 0;
}
