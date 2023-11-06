#include "../libremia.h"


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



void RandomWalk(bool cont,const unsigned int& K, vector<double>& dist, Random & rand){
	
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