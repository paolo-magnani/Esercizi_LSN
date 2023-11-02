#include "genalg.h"

int main(){ 

	/*Pop_circle pop(100, 34);
	//pop.print_pop();
	//pop.print_city();

	//cout << endl << endl;
	
	cout << endl << "In Evolve, mutazioni + crossover:" << endl;
	cout << endl;

	pop.evolveL1(10000);
	pop.getresults();*/

	Pop_square pop2(100, 34);
	pop2.print_pop();
	//pop.print_city();

	cout << endl << endl;
	
	cout << endl << "In Evolve, mutazioni + crossover:" << endl;
	cout << endl;

	pop2.evolveL1(100000);
	pop2.getresults();


	return 0;
}
