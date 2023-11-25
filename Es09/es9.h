#include "../libremia.h"

bool comp(const vector<double>& pippo, const vector<double> & pluto){ // funzione di comparazione per il sort
		
			if( pippo.back() < pluto.back() ) return true;
			else return false;
}

struct city_in_square{ //posizione per ogni città
	vector<double> x;
	vector<double> y;
};

class Pop_circle {
	
	private:
	
		vector<vector<double>> genome; //il genoma è un vettore di cromosomi (ognuno contenente una sequenza di città)
		vector<double> city; // sequenza contenente tutte le distanze delle città dalla prima
		vector<double> best, best_mean, best_chromo; // migliore distanza ad ogni generazione, migliore media su metà del genoma, miglior percorso
		unsigned int generations; //numero di generazioni
		double p_cross=0.5; // probabilità di crossing
		double p_mut=0.1; // probabilità di mutazione
		Random rand;
		bool dist;
		
	public:
		
		//Default constructor
		Pop_circle(){}
		
		//constructor
		Pop_circle(unsigned int n_chr, unsigned int n_genes){

			iniz(rand); // inizializzo il generatore random
				
			vector<double> chromo(n_genes);
			city.push_back(0.); //fisso la prima città
			for(unsigned int i=0; i<n_genes; i++){
				chromo[i]=i;
				if(i!=0) city.push_back(rand.Rannyu(0.,2*M_PI)); //genero le altre città
			}
			
			sort(city.begin(),city.end()); //per fare i test, la sequenza corretta più breve tra le città è 1-2-3-4-5-6 ...
			
			for(unsigned int i=0; i<n_chr; i++){	//preparo i cromosomi con sequenze casuali, aggiungo la misura e check per vedere se sono giusti
				random_shuffle(chromo.begin()+1, chromo.end());
				genome.push_back(chromo);
                genome[i].push_back(eval1(chromo));
                check(genome[i]);
			}

		}
		
		//Destructor
		~Pop_circle(){}
		
		//functions	

        void set_probs(double p_c, double p_m){
			
			p_cross=p_c;
			p_mut=p_m;
		}
		
		int selection(){ //questa è una selezione truccata in cui vengono scelti di più i primi della lista
			return int(genome.size()*pow(rand.Rannyu(),7.));
		}		
		
		void gensort(){ //uso il sort della standard library ma con una legge di scelta costruita ad hoc
			sort(genome.begin(), genome.end(), comp); 
		}
		
        bool check(const vector<double>& chr){ //sommati tutti insieme, i valori nei cromosomi (dotati di misura) devono dare N*(N-1)/2 
			double sum = 0;
			for(unsigned int i=0; i<chr.size()-1; i++) sum += chr[i];
			double realvalue = (city.size()-1)*city.size()/2;
			if(realvalue!=sum){ 
				cerr << endl << "WARNING: chromosome is not correct" << endl;
				for(unsigned int i=0;  i<chr.size(); i++) cout << chr[i] << " ";
				cout << endl;
				return true;
			}
			else return false;
		}
		
		double eval1(const vector<double>& chromo){ //valuto la distanza per per un vettore di città posizionate su un cerchio
			double L=0., Lx=0., Ly=0.;
			unsigned int index1;
			unsigned int index2;
			
			for(unsigned int i=0; i<(chromo.size()-1); i++){
				index1 = chromo[i];
				index2 = chromo[i+1];
				Lx = cos(city[index1])-cos(city[index2]);
				Ly = sin(city[index1])-sin(city[index2]);
				L += sqrt((Lx*Lx)+(Ly*Ly));
			}
			unsigned int appo = chromo.size()-1;
			index2 = chromo[appo];
			
			Lx=cos(city[0])-cos(city[index2]);
			Ly=sin(city[0])-sin(city[index2]);
			L += sqrt((Lx*Lx)+(Ly*Ly));
			
			return L;
		}
		
        double distance(const vector<double>& chromo){ //valuto la distanza sui cromosomi
			double L=0., Lx=0., Ly=0.;
			unsigned int index1;
			unsigned int index2;
			
			for(unsigned int i=0; i<(chromo.size()-2); i++){
				index1 = chromo[i];
				index2 = chromo[i+1];
				Lx = cos(city[index1])-cos(city[index2]);
				Ly = sin(city[index1])-sin(city[index2]);
				L += sqrt((Lx*Lx)+(Ly*Ly));
			}
			unsigned int appo = chromo.size()-2;
			index2 = chromo[appo];
			
			Lx=cos(city[0])-cos(city[index2]);
			Ly=sin(city[0])-sin(city[index2]);
			L += sqrt((Lx*Lx)+(Ly*Ly));
			
			return L;            
        }

		void setL1(){ //setto la distanza del tragitto per ciascun cromosoma (se non è stata già settata)
            for(unsigned int i=0; i<genome.size(); i++){
				if(genome[i].size()!=city.size()+1) genome[i].push_back(eval1(genome[i]));
			}
		}
		
        void print_pop(){ //stampo tutta il genoma composto dai cromosomi
            cout << "Stampo tutto il genoma: " << endl << endl;
			for(unsigned int j=0; j<genome.size(); j++){
                cout << "Cromosoma n° " << j + 1 << ": " ;
                for(unsigned int i=0; i<genome[j].size(); i++) cout << genome[j][i] << " ";
                cout << endl << endl;
			}

		}

        void print_city(){ //stampo il vettore di città
			for(unsigned int i=0; i<city.size(); i++) cout << i << " :  " << city[i] << endl;
		}
        
        void evolveL1(unsigned int n_rep){ //evolvo la popolazione e registro i migliori
			
			generations = n_rep;
			gensort();
            print_pop();
            for(unsigned int i = 0; i<n_rep; i++){
                crossover();
				mutations();
                gensort();

				best.push_back(genome[0].back());

				double media=0;
				for(unsigned int j=0; j<genome.size()/2; j++) media += genome[j].back();
				best_mean.push_back(2.*media/double(genome.size()));

				if((i+1)%int(n_rep/20)==0) cout << endl << "Progresso: " << (i+1)*100/n_rep << "%   Distanza migliore: " << genome[0].back();
            }
			cout << endl << endl;
			best_chromo = genome[0];
            print_pop();		
		}

        void mutations(){ //contiene le mutazioni che subisce il genoma
			
			for(unsigned int chr=0; chr<genome.size();chr++){
			
				if(rand.Rannyu()<p_mut){//scambio due città 
					
					unsigned int idx, idx2;
					idx = rand.Rannyu(1., genome[chr].size()-1);
					
					do{ 
						idx2 = rand.Rannyu(1., genome[chr].size()-1); 
					}while(idx==idx2);
					
					swap(genome[chr][idx], genome[chr][idx2]);
					genome[chr].back() = distance(genome[chr]);
					bool sentinel = check(genome[chr]);
					if(sentinel) cout << " Failure in mutation 1 " << endl << endl;
				}	
				

				if(rand.Rannyu()<p_mut){//sposto un blocco di n città di m posizioni 
					unsigned int idx = rand.Rannyu(1., genome[chr].size()-2);
					//cout << endl << "idx = " << idx << endl; 
					unsigned int m = rand.Rannyu(1., genome[chr].size()-1-idx); 
					//cout << endl << "m = " << m << endl; 
					unsigned int n = rand.Rannyu(1.,genome[chr].size()-1-idx-m);
					//cout << endl << "n = " << n << endl;
					vector<double> appo;
					for(unsigned int i=0; i<n; i++) appo.push_back(genome[chr][idx+i]);
					for(unsigned int i=0; i<m; i++)	genome[chr][idx+i] = genome[chr][idx+n+i];
					for(unsigned int i=0; i<n; i++) genome[chr][idx+m+i] = appo[i];
					genome[chr].back() = distance(genome[chr]);
					bool sentinel = check(genome[chr]);
					if(sentinel) cout << " Failure in mutation 2 " << endl << endl;
				}
				
				if(rand.Rannyu()<p_mut){//inverto l'ordine di n città
					unsigned int idx = rand.Rannyu(1., genome[chr].size()-2);
					unsigned int n = rand.Rannyu(1., genome[chr].size()-1-idx);
					vector<double> appo;
					for(unsigned int i=0; i<n; i++) appo.push_back(genome[chr][idx+i]);	
					for(unsigned int i=0; i<n; i++) genome[chr][idx+i] = appo[n-1-i];
					genome[chr].back() = distance(genome[chr]);
					bool sentinel = check(genome[chr]);
					if(sentinel) cout << " Failure in mutation 3 " << endl << endl;
				}
				
				
				if(rand.Rannyu()<p_mut){//permuto l'ordine di n città
					int idx = rand.Rannyu(1., genome[chr].size()-2);
					int n = rand.Rannyu(1., genome[chr].size()-1-idx);
					random_shuffle(&genome[chr][idx], &genome[chr][idx+n]);
					genome[chr].back() = distance(genome[chr]);
					bool sentinel = check(genome[chr]);
					if(sentinel) cout << " Failure in mutation 4 " << endl << endl;
				}

			}

		}
		
		void crossover(){ //strong crossover
			vector<vector<double>> prole;

			for(unsigned int k=0; k<genome.size()/2; k++){

				int dad = selection();
				int mum = selection();
					
				if(rand.Rannyu()<p_cross){
				
					unsigned int idx = rand.Rannyu(1., genome[mum].size()-1); //è un matriarcato
					
					vector<double> dad_tail, mum_tail, son, daughter;
					
					for(unsigned int i=1; i<genome[mum].size()-1; i++){
						bool m_check = false, d_check = false;
						for(unsigned int j=1; j<idx; j++){
							if(genome[dad][i]==genome[mum][j]) m_check = true;
							if(genome[mum][i]==genome[dad][j]) d_check = true;
						}
						if(m_check == false) mum_tail.push_back(genome[dad][i]);
						if(d_check == false) dad_tail.push_back(genome[mum][i]);
					}
					
					if(mum_tail.size()!=genome[mum].size()-idx-1) cout << " there is a mum problem";
					if(dad_tail.size()!=genome[mum].size()-idx-1) cout << " there is a dad problem";
					
					
					auto mumcut = genome[mum].begin() + idx;
					auto dadcut = genome[dad].begin() + idx;

					// prima parte
					daughter.insert(daughter.begin(), genome[mum].begin(), mumcut);
					son.insert(son.begin(), genome[dad].begin(), dadcut);

					// seconda parte
					daughter.insert(daughter.end(), mum_tail.begin(), mum_tail.end());
					son.insert(son.end(), dad_tail.begin(), dad_tail.end());

					daughter.push_back(eval1(daughter));
					son.push_back(eval1(son));

					bool sent1 = check(daughter);
					bool sent2 = check(son);
					if(sent1){
						cout << " Failure in crossover: mum will fix it " << endl << endl;
						daughter=genome[mum];
					}
					if(sent2){
						cout << " Failure in crossover: dad will fix it " << endl << endl;
						son=genome[dad];
					}
					prole.push_back(daughter);
					prole.push_back(son);
					
				}
				else{
					
					prole.push_back(genome[mum]);
					prole.push_back(genome[dad]);
				}
			}

			genome = prole;
			gensort();
		}

		void getresults(){ //stampo risultati su file
			ofstream outbest, outbest_chromo, out_evolution;
			outbest.open("best_distance.dat");
			outbest_chromo.open("best_path.dat");
			out_evolution.open("risultati/best_evolution.dat");

			for(unsigned int i = 0; i < generations; i++) outbest << i+1 << "," << best[i] << "," << best_mean[i] << endl;
			for(unsigned int i = 0; i < best_chromo.size()-1; i++){
				double x = cos(city[best_chromo[i]]);
				double y = sin(city[best_chromo[i]]);
				outbest_chromo << i << "," << x << "," << y << endl;
			}
			outbest.close();
			outbest_chromo.close();
			out_evolution.close();
		}		
		
};

class Pop_square {
	
	private:
	
		vector<vector<double>> genome; //il genoma è un vettore di cromosomi (ognuno contenente una sequenza di città)
		city_in_square city; // sequenza contenente tutte le distanze delle città dalla prima
		vector<double> best, best_mean, best_chromo;
		unsigned int generations; //numero di generazioni
		double p_cross=0.5; // probabilità di crossing
		double p_mut=0.1; // probabilità di mutazione
		Random rand;
		bool dist;
		
	public:
		
		//Default constructor
		Pop_square(){}
		
		//constructor
		Pop_square(unsigned int n_chr, unsigned int n_genes){

			iniz(rand);
				
			vector<double> chromo(n_genes);
			city.x.push_back(0.);
			city.y.push_back(0.); //fisso le coordinate della prima città
			for(unsigned int i=0; i<n_genes; i++){
				chromo[i]=i;
				if(i!=0){
					city.x.push_back(rand.Rannyu(0.,1.));
					city.y.push_back(rand.Rannyu(0.,1.)); //genero le altre città
				}
			}
			
			for(unsigned int i=0; i<n_chr; i++){	//preparo i cromosomi con sequenze casuali, aggiungo la misura e check per vedere se sono giusti
				random_shuffle(chromo.begin()+1, chromo.end());
				genome.push_back(chromo);
                genome[i].push_back(eval1(chromo));
                check(genome[i]);
			}

		}

		Pop_square(unsigned int n_chr, string filename){ // costruttore nel caso le città siano su un file

			iniz(rand);
			
			ifstream in;
			in.open(filename);
			double appo, appox, appoy;
			
			while(in >> appo >> appox >> appoy){
				city.x.push_back(appox);
				city.y.push_back(appoy);
			}

			vector<double> chromo(city.x.size());

			for(unsigned int i=0; i<city.x.size(); i++){
				chromo[i]=i;
			}
			
			for(unsigned int i=0; i<n_chr; i++){	//preparo i cromosomi con sequenze casuali, aggiungo la misura e check per vedere se sono giusti
				random_shuffle(chromo.begin()+1, chromo.end());
				genome.push_back(chromo);
                genome[i].push_back(eval1(chromo));
                check(genome[i]);
			}
			in.close();
		}
		
		//Destructor
		~Pop_square(){}
		
		//functions	

        void set_probs(double p_c, double p_m){
			
			p_cross=p_c;
			p_mut=p_m;
		}
		
		int selection(){ //questa è una selezione truccata in cui vengono scelti di più i primi della lista
			return int(genome.size()*pow(rand.Rannyu(),7.));
		}		
		
		void gensort(){ //uso il sort della standard library ma con una legge di scelta costruita ad hoc
			sort(genome.begin(), genome.end(), comp); 
		}
		
        bool check(const vector<double>& chr){ //sommati tutti insieme, i valori nei cromosomi (dotati di misura) devono dare N*(N-1)/2 
			double sum = 0;
			for(unsigned int i=0; i<chr.size()-1; i++) sum += chr[i];
			double realvalue = (city.x.size()-1)*city.x.size()/2;
			if(realvalue!=sum){ 
				cerr << endl << "WARNING: chromosome is not correct" << endl;
				for(unsigned int i=0;  i<chr.size(); i++) cout << chr[i] << " ";
				cout << endl;
				return true;
			}
			else return false;
		}
		
		double eval1(const vector<double>& chromo){ //valuto la distanza per per un vettore di città posizionate in quadrato
			double L=0., Lx=0., Ly=0.;
			unsigned int index1;
			unsigned int index2;
			
			for(unsigned int i=0; i<(chromo.size()-1); i++){
				index1 = chromo[i];
				index2 = chromo[i+1];
				Lx = city.x[index1]-city.x[index2];
				Ly = city.y[index1]-city.y[index2];
				L += sqrt((Lx*Lx)+(Ly*Ly));
			}
			unsigned int appo = chromo.size()-1;
			index2 = chromo[appo];
			
			Lx=city.x[0]-city.x[index2];
			Ly=city.y[0]-city.y[index2];
			L += sqrt((Lx*Lx)+(Ly*Ly));
			
			return L;
		}
		
        double distance(const vector<double>& chromo){ //valuto la distanza sui cromosomi
			double L=0., Lx=0., Ly=0.;
			unsigned int index1;
			unsigned int index2;
			
			for(unsigned int i=0; i<(chromo.size()-2); i++){
				index1 = chromo[i];
				index2 = chromo[i+1];
				Lx = city.x[index1]-city.x[index2];
				Ly = city.y[index1]-city.y[index2];
				L += sqrt((Lx*Lx)+(Ly*Ly));
			}
			unsigned int appo = chromo.size()-2;
			index2 = chromo[appo];
			
			Lx=city.x[0]-city.x[index2];
			Ly=city.y[0]-city.y[index2];
			L += sqrt((Lx*Lx)+(Ly*Ly));
			
			return L;            
        }

		void setL1(){ //setto la distanza del tragitto per ciascun cromosoma (se non è stata già settata)
            for(unsigned int i=0; i<genome.size(); i++){
				if(genome[i].size()!=city.x.size()+1) genome[i].push_back(eval1(genome[i]));
			}
		}
		
        void print_pop(){ //stampo tutta il genoma composto dai cromosomi
            cout << "Stampo tutto il genoma: " << endl << endl;
			for(unsigned int j=0; j<genome.size(); j++){
                cout << "Cromosoma n° " << j + 1 << ": " ;
                for(unsigned int i=0; i<genome[j].size(); i++) cout << genome[j][i] << " ";
                cout << endl << endl;
			}

		}

        void print_city(){ //stampo il vettore di città
			for(unsigned int i=0; i<city.x.size(); i++) cout << "città " << i << " :  " << city.x[i] << " , " << city.y[i] << endl;
		}
        
        void evolveL1(unsigned int n_rep){ //evolvo la popolazione e registro i migliori
			
			generations = n_rep;
			gensort();
            print_pop();
            for(unsigned int i = 0; i<n_rep; i++){
                crossover();
				gensort();
				mutations();
                gensort();

				best.push_back(genome[0].back());

				double media=0;
				for(unsigned int j=0; j<genome.size()/2; j++) media += genome[j].back();
				best_mean.push_back(2.*media/double(genome.size()));

				if((i+1)%int(n_rep/20)==0) cout << endl << "Progresso: " << (i+1)*100/n_rep << "%   Distanza migliore: " << genome[0].back();
            }
			cout << endl << endl;
			best_chromo = genome[0];
            print_pop();		
		}

        void mutations(){ //contiene le mutazioni che subisce il genoma
			
			for(unsigned int chr=0; chr<genome.size();chr++){
			
				if(rand.Rannyu()<p_mut){//scambio due città 
					
					unsigned int idx, idx2;
					idx = rand.Rannyu(1., genome[chr].size()-1);
					
					do{ 
						idx2 = rand.Rannyu(1., genome[chr].size()-1); 
					}while(idx==idx2);
					
					swap(genome[chr][idx], genome[chr][idx2]);
					genome[chr].back() = distance(genome[chr]);
					bool sentinel = check(genome[chr]);
					if(sentinel) cout << " Failure in mutation 1 " << endl << endl;
				}	


				if(rand.Rannyu()<p_mut){//sposto un blocco di n città di m posizioni 
					unsigned int idx = rand.Rannyu(1., genome[chr].size()-2);
					unsigned int m = rand.Rannyu(1., genome[chr].size()-1-idx); 
					unsigned int n = rand.Rannyu(1.,genome[chr].size()-1-idx-m);

					vector<double> appo;
					for(unsigned int i=0; i<n; i++) appo.push_back(genome[chr][idx+i]);
					for(unsigned int i=0; i<m; i++)	genome[chr][idx+i] = genome[chr][idx+n+i];
					for(unsigned int i=0; i<n; i++) genome[chr][idx+m+i] = appo[i];
					genome[chr].back() = distance(genome[chr]);
					bool sentinel = check(genome[chr]);
					if(sentinel) cout << " Failure in mutation 2 " << endl << endl;
				}

				
				if(rand.Rannyu()<p_mut){//inverto l'ordine di n città
					unsigned int idx = rand.Rannyu(1., genome[chr].size()-2);
					unsigned int n = rand.Rannyu(1., genome[chr].size()-1-idx);
					vector<double> appo;
					for(unsigned int i=0; i<n; i++) appo.push_back(genome[chr][idx+i]);	
					for(unsigned int i=0; i<n; i++) genome[chr][idx+i] = appo[n-1-i];
					genome[chr].back() = distance(genome[chr]);
					bool sentinel = check(genome[chr]);
					if(sentinel) cout << " Failure in mutation 3 " << endl << endl;
				}

				
				if(rand.Rannyu()<p_mut){//permuto l'ordine di n città
					int idx = rand.Rannyu(1., genome[chr].size()-2);
					int n = rand.Rannyu(1., genome[chr].size()-1-idx);
					random_shuffle(&genome[chr][idx], &genome[chr][idx+n]);
					genome[chr].back() = distance(genome[chr]);
					bool sentinel = check(genome[chr]);
					if(sentinel) cout << " Failure in mutation 4 " << endl << endl;
				}
			}
		}
		
		void crossover(){ //strong crossover
			vector<vector<double>> prole;

			for(unsigned int k=0; k<genome.size()/2; k++){

				int dad = selection();
				int mum = selection();
					
				if(rand.Rannyu()<p_cross){
				
					unsigned int idx = rand.Rannyu(1., genome[mum].size()-1); //è un matriarcato
					
					vector<double> dad_tail, mum_tail, son, daughter;
					
					for(unsigned int i=1; i<genome[mum].size()-1; i++){
						bool m_check = false, d_check = false;
						for(unsigned int j=1; j<idx; j++){
							if(genome[dad][i]==genome[mum][j]) m_check = true;
							if(genome[mum][i]==genome[dad][j]) d_check = true;
						}
						if(m_check == false) mum_tail.push_back(genome[dad][i]);
						if(d_check == false) dad_tail.push_back(genome[mum][i]);
					}
					
					if(mum_tail.size()!=genome[mum].size()-idx-1) cout << " there is a mum problem";
					if(dad_tail.size()!=genome[mum].size()-idx-1) cout << " there is a dad problem";
					
					
					auto mumcut = genome[mum].begin() + idx;
					auto dadcut = genome[dad].begin() + idx;

					// prima parte
					daughter.insert(daughter.begin(), genome[mum].begin(), mumcut);
					son.insert(son.begin(), genome[dad].begin(), dadcut);

					// seconda parte
					daughter.insert(daughter.end(), mum_tail.begin(), mum_tail.end());
					son.insert(son.end(), dad_tail.begin(), dad_tail.end());

					daughter.push_back(eval1(daughter));
					son.push_back(eval1(son));

					bool sent1 = check(daughter);
					bool sent2 = check(son);
					if(sent1){
						cout << " Failure in crossover: mum will fix it " << endl << endl;
						daughter=genome[mum];
					}
					if(sent2){
						cout << " Failure in crossover: dad will fix it " << endl << endl;
						son=genome[dad];
					}
					prole.push_back(daughter);
					prole.push_back(son);
					
				}
				else{
					
					prole.push_back(genome[mum]);
					prole.push_back(genome[dad]);
				}
			}

			genome = prole;
			gensort();
		}

		void getresults(){ //stampo risultati su file
			ofstream outbest, outbest_chromo;
			outbest.open("best_distance_square.dat");
			outbest_chromo.open("best_path_square.dat");

			for(unsigned int i = 0; i < generations; i++) outbest << i+1 << "," << best[i] << "," << best_mean[i] << endl;
			for(unsigned int i = 0; i < best_chromo.size()-1; i++){
				double x = city.x[best_chromo[i]];
				double y = city.y[best_chromo[i]];
				outbest_chromo << i << "," << x << "," << y << endl;
			}
			outbest.close();
			outbest_chromo.close();
		}		
		
};
