#include "es9.h"

int main(){ 

	//miglior percorso per città su una circonferenza
	Pop_circle pop(100, 34); // 100 cromosomi, 34 città

	cout << endl << endl << " Cerco il migliore percorso per le città posizionate sulla circonferenza... " << endl << endl;
	pop.evolveL1(10000);
	pop.getresults();


	//miglior percorso per città in un quadrato
	cout << endl << endl << " Città in un quadrato: " << endl;
	Pop_square pop2(100, 34); // 100 cromosomi, 34 città

	cout << endl << endl << " Cerco il migliore percorso per le città posizionate in un quadrato... " << endl << endl;
	pop2.evolveL1(10000);
	pop2.getresults();

	return 0;
}
