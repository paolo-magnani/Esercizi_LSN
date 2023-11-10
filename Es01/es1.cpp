#include "../libremia.h" //le funzioni o classi utilizzate in più di una esercitazione sono state inserite in una libreria
#define M 1000000 // numero di lanci che intendo fare
#define N 100 // numero di step
#define nchi 10000
#define Mchi 100
#define lambda 1 //valori da inserire per le distribuzioni esponenziali e lorentziane
#define mu 0
#define gamma 1


int main(){
	
	const unsigned int L=M/N;
	Random rand; 
	iniz(rand); // inizializzo il generatore di numeri casuali
	
	vector<double> r(M), x(N);
	
	generate(r.begin(),r.end(),[&](){ return rand.Rannyu();}); // riempio il vettore 'r' di numeri random
	for(unsigned int i=0; i<N; i++) x[i]=i*L;
	
	//punto 1.1
	cout << " Risolvo il punto 1.1... " << endl;
	
	vector<double> su1_prog(N,0), er1_prog(N,0); // inizializzo i vector a 0
	
	mediablocchi(M, N, r, su1_prog, er1_prog); // calcolo la media a blocchi

	stampamedia("risultati/1011.dat", M, su1_prog, er1_prog);
	
	
	//punto 1.2
	//analogo del punto 1
	cout << " Risolvo il punto 1.2... " << endl;

	vector<double> su2_prog(N,0), er2_prog(N,0);
	
	for(unsigned int i=0; i<M; i++) r[i]= pow(r[i]-0.5, 2.);
	
	mediablocchi(M, N, r, su2_prog, er2_prog);
			
	stampamedia("risultati/1012.dat", M, su2_prog, er2_prog);
	
	//punto 1.3
	ofstream out;
	out.open("risultati/1013.dat");
	iniz(rand);	//riinizializzo il generatore random
	
	vector<double> chisq(Mchi,0); //vettore di risultati del chi quadro
	
	double low = 0., high=1.;  // definisco il valore minimo e massimo dell'intervallo
	
	for(unsigned int i=0; i<Mchi; i++){ // ripeto il processo Mchi volte
		
		vector<double> ran(nchi);
		generate(ran.begin(),ran.end(),[&](){ return rand.Rannyu();}); 	// riempio il vettore di numeri casuali estratti tra [0,1)
		
		vector<unsigned int> n = binsearch(ran, low, high, Mchi); // algoritmo di ricerca binaria che cicla su tutti i numeri generati per incasellarli nel bin corrispondente
		
		for(unsigned int j=0; j<Mchi; j++){
		
			chisq[i] += pow(double(n[j])-double(nchi/Mchi), 2.); // calcolo il chi quadro
		} 
		chisq[i] /= double(nchi/Mchi);
		out << chisq[i] << endl; // stampo su file
	}
	
	out.close();
	
	
	
	//punto 2: test teorema limite centrale
	cout << " Risolvo il punto 2..." << endl;
	ofstream outUni;
	ofstream outExp;
	ofstream outLor;
	outUni.open("risultati/Uni.dat");
	outExp.open("risultati/Exp.dat");
	outLor.open("risultati/Lor.dat");
	iniz(rand); 
	vector<double> uni(nchi,0), exp(nchi,0), lor(nchi,0);
	vector<double> rep = {1, 2, 10, 100}; // ripetizioni per il test del teorema del limite centrale
	
	for(unsigned int i=0; i<nchi; i++){ // riempio i vettori di 'nchi' occorrenze della media 
	
		int conta = 0;
		
		for(unsigned int j=0; j<4; j++){
		
			while(conta < rep[j]){ // sommo numeri estratti secondo le distribuzioni
				uni[i] += rand.Rannyu();
				exp[i] += rand.Exp(lambda);
				lor[i] += rand.Lorentz(mu, gamma);
				conta ++;
			}
			
			outUni << uni[i]/rep[j] << "," ; // calcolo la media
			outExp << exp[i]/rep[j] << "," ;
			outLor << lor[i]/rep[j] << "," ;
		
		}
		
		outUni << endl;
		outExp << endl;
		outLor << endl;
	
	}
	
	outUni.close();
	outExp.close();
	outLor.close();
	
	
	
	//punto 3: ago di Buffon
	cout << " Risolvo il punto 3... " << endl;
	iniz(rand); //riinizializzo il generatore di numeri random 
		
	const double d = 1.; // dimensione del box
	const double lenght = 0.4; // lunghezza dell'ago di Buffon
	unsigned int Nhit = 0; // numero di intersezioni tra l'ago e i confini del box
	
	vector<double> su3(N,0), su3sq(N,0), er3(N,0);
	double sum = 0, sumsq=0;
	
	
	for(unsigned int i=0; i<N; i++){
	
		for(unsigned int j=0; j<L; j++){ 
		
			double x0 = rand.Rannyu(0.,d); // l'ago cade generica posizione x

			double x1 = rand.Rannyu(-1.,1.); // uso metodo Accept-Reject per campionare l'angolo
			double y1 = rand.Rannyu(-1.,1.);
			
			while(sqrt(pow(x1,2.)+pow(y1,2.))>1.){ // prendo punti all'interno di un cerchio unitario
				x1 = rand.Rannyu(-1.,1.);
				y1 = rand.Rannyu(-1.,1.);	
			}
			
			double sin = x1/sqrt(pow(x1,2.)+pow(y1,2.));
			
			if(x0+(sin*lenght)<0 or x0+(sin*lenght)>d) Nhit++; // se l'ago interseca i confini del box, aumenta 'Nhit'
		}
		
		sum += (2.*lenght*L)/(Nhit*d); // stima di Buffon
		su3[i] = sum/(i+1.);
		
		sumsq += pow((2.*lenght*L)/(Nhit*d),2.);
		su3sq[i] = sumsq/(i+1.);
		
		er3[i]=error(su3, su3sq, i);
		
		Nhit=0;
	}

	stampamedia("risultati/Buffon.dat", M, su3, er3);
	
	cout << " Completo." << endl;
	 
	 
	
 	return 0;
}
