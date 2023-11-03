#include "../libremia.h" //le funzioni o classi utilizzate in più di una esercitazione sono state inserite in una cartella superiore
#define M 1000000 // definisco il numero di lanci che intendo fare
#define N 100 // definisco il numero di step
#define nchi 10000
#define Mchi 100
#define lambda 1 //valori da inserire per le ditribuzioni esponenziali e lorentziane
#define mu 0
#define gamma 1


int main(){
	
	const unsigned int L=M/N;
	Random rand; 
	iniz(rand); // inizializzo il generatore di numeri casuali
	
	vector<double> r(M), x(N);
	
	generate(r.begin(),r.end(),[&](){ return rand.Rannyu();}); // uso una lambda function per riempire il vettore 'r' di numeri generati randomicamente
	for(unsigned int i=0; i<N; i++) x[i]=i*L;
	
	//punto 1.1
	cout << " Risolvo il punto 1.1... " << endl;
	
	vector<double> su1_prog(N,0), er1_prog(N,0); // uso un costruttore della classe vector in modo da creare un vettore di zeri
	
	mediablocchi(M, N, r, su1_prog, er1_prog); // uso la funzione contenuta nella libreria

	stampamedia("1011.dat", M, su1_prog, er1_prog);
	
	
	//punto 1.2
	//analogo del punto 1
	cout << " Risolvo il punto 1.2... " << endl;

	vector<double> su2_prog(N,0), er2_prog(N,0);
	
	for(unsigned int i=0; i<M; i++) r[i]= pow(r[i]-0.5, 2.);
	
	mediablocchi(M, N, r, su2_prog, er2_prog);
			
	stampamedia("1012.dat", M, su2_prog, er2_prog);
	
	//punto 1.3
	ofstream out;
	out.open("1013.dat");
	iniz(rand);	//riinizializzo il generatore random
	
	vector<double> chisq(Mchi,0); //vettore di risultati del chi quadro inizializzato a 0
	
	double ampl = 1./double(Mchi); // definisco l'ampiezza di ogni bin
	double low = 0.;  // definisco il valore minimo dell'intervallo
	
	
	for(unsigned int i=0; i<Mchi; i++){ // ripeto il processo Mchi volte
		
		vector<double> ran(nchi);
		generate(ran.begin(),ran.end(),[&](){ return rand.Rannyu();}); 	// riempio il vettore di numeri casuali estratti tra [0,1)
		vector<int> n(Mchi,0); 	// creo il vettore che conta le occorrenze in ogni bin
		
		for(unsigned int j=0; j<nchi; j++){  // algoritmo di ricerca binaria che cicla su tutti i numeri generati per incasellarli nel bin corrispondente
			
			unsigned int n_min = 0; 
			unsigned int n_max = Mchi - 1;
		
			while(n_min<=n_max){
				
				unsigned int medio = (n_min + n_max)/2;
				double bin_min = double(low) + double(medio)*ampl;
				double bin_max = double(low) + double(medio+1)*ampl;
				
				if(ran[j]>=bin_min and ran[j]<bin_max){
					n[medio]++;
					break;
				}
				else if(ran[j]<bin_min) n_max = medio-1;
				else n_min = medio+1;
				}
		}
		
		for(unsigned int j=0; j<Mchi; j++){
		
			chisq[i] += pow(double(n[j])-double(nchi/Mchi), 2.); // calcolo il chi quadro
		} 
		chisq[i] /= double(nchi/Mchi);
		out << chisq[i] << endl; // stampo su file
	}
	
	out.close();
	
	
	
	//punto 2: test teorema limite centrale
	
	ofstream outUni;
	ofstream outExp;
	ofstream outLor;
	outUni.open("Uni.dat");
	outExp.open("Exp.dat");
	outLor.open("Lor.dat");
	iniz(rand); 
	vector<double> uni(nchi,0), exp(nchi,0), lor(nchi,0);
	vector<double> rep = {1, 2, 10, 100};
	
	for(unsigned int i=0; i<10000; i++){
	
		int conta = 0;
		
		for(unsigned int j=0; j<4; j++){
		
			while(conta < rep[j]){
				uni[i] += rand.Rannyu();
				exp[i] += rand.Exp(lambda);
				lor[i] += rand.Lorentz(mu, gamma);
				conta ++;
			}
			
			outUni << uni[i]/rep[j] << "," ;
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

	iniz(rand);
	ofstream out4;
	out4.open("Buffon.dat");
		
	const double d = 1.;
	const double lenght = 0.4;
	unsigned int Nhit = 0;
	
	vector<double> su3(N,0), su3sq(N,0), er3(N,0);
	double sum = 0, sumsq=0;
	
	
	for(unsigned int i=0; i<N; i++){
	
		for(unsigned int j=0; j<L; j++){				// L è la dimensione di ogni blocco 
		
			double x0 = rand.Rannyu(0.,d);
			double x1 = rand.Rannyu(-1.,1.);
			double y1 = rand.Rannyu(-1.,1.);
			
			while(sqrt(pow(x1,2.)+pow(y1,2.))>1.){
				x1 = rand.Rannyu(-1.,1.);
				y1 = rand.Rannyu(-1.,1.);	
			}
			
			double sin = x1/sqrt(pow(x1,2.)+pow(y1,2.));
			
			if(x0+(sin*lenght)<0 or x0+(sin*lenght)>d) Nhit++;	
		}
		
		sum += (2.*lenght*L)/(Nhit*d);
		su3[i] = sum/(i+1.);
		
		sumsq += pow((2.*lenght*L)/(Nhit*d),2.);
		su3sq[i] = sumsq/(i+1.);
		
		er3[i]=error(su3, su3sq, i);
		
		Nhit=0;
	}
	
	//mediablocchi(M, N, stima, su3, er3);
	
	for(unsigned int i=0; i<N; i++) out4 << i+1 << "," << su3[i] << "," << er3[i] << endl;
	
	out4.close();
	
	cout << endl << "Completo." << endl;
	 
	 
	
 	return 0;
}
