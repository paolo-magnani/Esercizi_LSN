#pragma once

#include "random/random.h"
#include <stdlib.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iomanip>

using namespace std;

double error(const vector<double> & AV,const vector<double> & AV2, unsigned int n){
	if(n==0) return 0;
	else return sqrt((AV2[n] - pow(AV[n],2.))/double(n));
}


void iniz(Random &rnd){

   int seed[4];
   int p1, p2;
   ifstream inprimes("../random/Primes");
   if (inprimes.is_open()){
      inprimes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   inprimes.close();

   ifstream input("../random/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   
}

/*void mediablocchi( const unsigned int M, const unsigned int N, const vector<double>& r, vector<double>& sum_prog, vector<double>& err_prog){

	vector<double> ave(N,0), av2(N,0), su2_prog(N,0);
	
	const unsigned int L=M/N;
	for(unsigned int i=0; i<N; i++){
		double sum = 0;
		
	for(unsigned int j=0; j<L; j++){
			double k=j+(i*L);
			sum += r[k];
			}
			
		ave[i] = sum/double(L);
		av2[i] = pow(ave[i], 2.);
	}
			
	for(unsigned int i=0; i<N; i++){
		for(unsigned int j=0; j<i+1; j++){
				sum_prog[i] += ave[j];
				su2_prog[i] += av2[j];
		}
		sum_prog[i] /= double(i+1.);
		su2_prog[i] /= double(i+1.);
		err_prog[i] = error(sum_prog, su2_prog, i);
	}

}*/ 

void mediablocchi( const unsigned int M, const unsigned int N, const vector<double>& r, vector<double>& sum_prog, vector<double>& err_prog){

	vector<double> su2_prog(N,0);
	double ave=0, av2=0;
	
	const unsigned int L=M/N;
	
	for(unsigned int i=0; i<N; i++){
		double sum = 0;
		
	for(unsigned int j=0; j<L; j++){
			double k=j+(i*L);
			sum += r[k];
			}
			
		ave += sum/double(L);
		av2 += pow(sum/double(L), 2.);
		
		sum_prog[i] = ave/double(i+1.);
		su2_prog[i] =av2/double(i+1.);
		err_prog[i] = error(sum_prog, su2_prog, i);
	}
			

} 
