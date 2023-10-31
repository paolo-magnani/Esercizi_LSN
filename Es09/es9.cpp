#include "genalg.h"

int main(){

	Pop pop(10, 34);
	pop.print_pop();
	//pop.print_city();

	cout << endl << endl;
	
	cout << endl << "In Evolve:" << endl;
	cout << endl;
	
	pop.evolveL1();
	//pop.evolveL1();
	//pop.gensort();

	return 0;
}
