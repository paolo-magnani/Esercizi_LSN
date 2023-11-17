//#include "es9.h"
#include "es9pointers.h"

int main(){ 

	// vector<double> z(10);
	// vector<vector<double>> v, w, x;
	// v.push_back(z);
	// x = v;
	// v.push_back(z);
	// auto it = v[0].begin()+9;
	// v[0][1]=2;
	// v[0][9]=4;
	// cout << endl << v.data()+1;
	// cout << endl << &v[1] << endl << endl;
	// cout << *(v[0].data()+1) << endl << endl;
	// cout << &it << endl << endl;
	// cout << &v[0].back() << endl << endl;

	// cout << endl << "Test con w"<< endl; 
	// w=v;
	// cout << endl << w.data()+1;
	// cout << endl << &w[1] << endl << endl;
	// cout << *(w[0].data()+1) << endl << endl;
	// cout << *(w[0].end()-1) << endl << endl;
	// cout << w[0].back() << endl << endl;
	// w=x;
	// cout << endl << w.data()+1;
	// cout << endl << &w[1] << endl << endl;
	// cout << *(w[0].data()+1) << endl << endl;
	// cout << *(w[0].end()-1) << endl << endl;
	// cout << w[0].back() << endl << endl;





	// Pop_circle pop(100, 34);
	// //pop.print_pop();
	// //pop.print_city();

	// //cout << endl << endl;
	
	// cout << endl << "In Evolve, mutazioni + crossover:" << endl;
	// cout << endl;

	// pop.evolveL1(10000);
	// pop.getresults();

	Pop_square pop2(100, "../../square.out");
	pop2.print_pop();
	//pop.print_city();

	cout << endl << endl;
	
	cout << endl << "In Evolve, mutazioni + crossover:" << endl;
	cout << endl;

	pop2.evolveL1(100000);
	pop2.getresults();


	return 0;
}
