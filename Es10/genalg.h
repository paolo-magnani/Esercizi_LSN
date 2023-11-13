#include "../libremia.h"
#include "mpi.h"

bool comp(vector<double> pippo, vector<double> pluto){
		
			if( pippo[pippo.size()-1]<pluto[pluto.size()-1]) return true;
			else return false;
}

struct american_city{ //posizione per ogni città
	//vector<string> name;
	//vector<string> state;
	vector<double> x;
	vector<double> y;
};

class Pop{
	
	private:
	
		vector<vector<double>> genome; //il genoma è un vettore di cromosomi (ognuno contenente una sequenza di città)
		american_city city; // sequenza contenente tutte le distanze delle città dalla prima
		vector<double> best, best_mean, best_chromo;
		unsigned int generations; //numero di generazioni
		double p_cross=0.5; // probabilità di crossing
		double p_mut=0.2; // probabilità di mutazione
		double p_mut2=0.3; // probabilità di scambio tra due città
		Random rand;
		bool dist;
		
	public:
		
		//Default constructor
		Pop(){}
		
		//constructor
		Pop(unsigned int n_chr, american_city capitals){

			iniz(rand);
				
			vector<double> chromo(capitals.x.size());
			
			for(unsigned int i=0; i<capitals.x.size(); i++){
				chromo[i]=i;
				city.x.push_back(capitals.x[i]);
				city.y.push_back(capitals.y[i]); // potrei anche uguagliare le due struct
			}
			
			for(unsigned int i=0; i<n_chr; i++){	//preparo i cromosomi con sequenze casuali, aggiungo la misura e check per vedere se sono giusti
				random_shuffle(chromo.begin()+1, chromo.end());
				genome.push_back(chromo);
                genome[i].push_back(eval1(chromo));
                check(genome[i]);
			}

		}
		
		//Destructor
		~Pop(){}
		
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
		
		double eval1(const vector<double>& chromo){ //valuto la distanza per per un vettore di città posizionate sul suolo americano
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
			unsigned int appo = chromo.size()-1;
			index2 = chromo[appo];
			
			//if(city[index2]>M_PI) L += fabs(2*M_PI-city[index2]);
			//else L += fabs(city[index2]);
			
			Lx=city.x[0]-city.x[index2];
			Ly=city.y[0]-city.y[index2];
			L += sqrt((Lx*Lx)+(Ly*Ly));
			
			return L;
		}
		
        double distance(const vector<double>& chromo){ //valuto la distanza sui cromosomi
			double L=0., Lx=0., Ly=0.;
			unsigned int index1;
			unsigned int index2;
			
			for(unsigned int i=0; i<(chromo.size()-3); i++){
				index1 = chromo[i];
				index2 = chromo[i+1];
				Lx = city.x[index1]-city.x[index2];
				Ly = city.y[index1]-city.y[index2];
				L += sqrt((Lx*Lx)+(Ly*Ly));
			}
			unsigned int appo = chromo.size()-2;
			index2 = chromo[appo];
			
			//if(city[index2]>M_PI) L += fabs(2*M_PI-city[index2]);
			//else L += fabs(city[index2]);
			
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
			for(unsigned int i=0; i<city.x.size(); i++) cout << i+1 << " :  " << city.x[i] << " , " << city.y[i] << endl;
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

				best.push_back(genome[0][genome[0].size()-1]);

				double media=0;
				for(unsigned int j=0; j<genome.size()/2; j++) media += genome[j][genome[j].size()-1];
				best_mean.push_back(2.*media/double(genome.size()));
            }
			cout << endl << endl;
			best_chromo = genome[0];
            print_pop();		
		}

        void mutations(){ //contiene le mutazioni che subisce il genoma
			
			int chr = selection();
			
			if(rand.Rannyu()<p_mut2){//scambio due città 
				
                int idx, idx2;
				idx = rand.Rannyu(1., genome[chr].size()-1);
				
                do{ 
                    idx2 = rand.Rannyu(1., genome[chr].size()-1); 
                }while(idx==idx2);
				
                swap(genome[chr][idx], genome[chr][idx2]);
				genome[chr][genome[chr].size()-1] = distance(genome[chr]);
				bool sentinel = check(genome[chr]);
				if(sentinel) cout << " Failure in mutation 1 " << endl << endl;
			}	
			
			chr = selection();

			if(rand.Rannyu()<p_mut){//sposto un blocco di n città di m posizioni 
				int idx = rand.Rannyu(1., genome[chr].size()-2);
				//cout << endl << "idx = " << idx << endl; 
				int m = rand.Rannyu(1., genome[chr].size()-1-idx); 
				//cout << endl << "m = " << m << endl; 
				int n = rand.Rannyu(1.,genome[chr].size()-1-idx-m);
				//cout << endl << "n = " << n << endl;
				vector<double> appo;
				for(unsigned int i=0; i<n; i++) appo.push_back(genome[chr][idx+i]);
				for(unsigned int i=0; i<m; i++)	genome[chr][idx+i] = genome[chr][idx+n+i];
				for(unsigned int i=0; i<n; i++) genome[chr][idx+m+i] = appo[i];
				genome[chr][genome[chr].size()-1] = distance(genome[chr]);
                bool sentinel = check(genome[chr]);
				if(sentinel) cout << " Failure in mutation 2 " << endl << endl;
			}

			chr = selection();
			
			if(rand.Rannyu()<p_mut){//inverto l'ordine di n città
				int idx = rand.Rannyu(1., genome[chr].size()-2);
				int n = rand.Rannyu(1., genome[chr].size()-1-idx);
				vector<double> appo;
				for(unsigned int i=0; i<n; i++) appo.push_back(genome[chr][idx+i]);	
				for(unsigned int i=0; i<n; i++) genome[chr][idx+i] = appo[n-1-i];
				genome[chr][genome[chr].size()-1] = distance(genome[chr]);
                bool sentinel = check(genome[chr]);
				if(sentinel) cout << " Failure in mutation 3 " << endl << endl;
			}
            
			chr = selection();
			
			if(rand.Rannyu()<p_mut){//permuto l'ordine di n città
				int idx = rand.Rannyu(1., genome[chr].size()-2);
				int n = rand.Rannyu(1., genome[chr].size()-1-idx);
				random_shuffle(&genome[chr][idx], &genome[chr][idx+n]);
				genome[chr][genome[chr].size()-1] = distance(genome[chr]);
                bool sentinel = check(genome[chr]);
				if(sentinel) cout << " Failure in mutation 4 " << endl << endl;
			}

		}
		
		void crossover(){ //accoppio le i migliori per generare prole ancora più bella
			
			if(rand.Rannyu()<p_cross){
			
				int dad = selection();
				int mum = selection();
				
				int idx = rand.Rannyu(1., genome[mum].size()-1); //è un matriarcato
				
				vector<double> dad_tail, mum_tail;
				
				for(unsigned int i=1; i<genome[mum].size()-1; i++){
					bool m_check = false, d_check = false;
					for(unsigned int j=1; j<idx; j++){
						if(genome[dad][i]==genome[mum][j]) m_check = true;
						if(genome[mum][i]==genome[dad][j]) d_check = true;
					}
					if(m_check == false) mum_tail.push_back(genome[dad][i]);
					if(d_check == false) dad_tail.push_back(genome[mum][i]);
				}
				
				if(mum_tail.size()!=genome[mum].size()-idx-1) cout << " UNO PROBLEMA!";
				if(dad_tail.size()!=genome[mum].size()-idx-1) cout << " DUO PROBLEMA!";
				
				
				if(mum_tail.size()==dad_tail.size()){ //rimpiazzo gli ultimi della lista
				
					for(unsigned int i=1; i<idx; i++){
						genome[genome.size()-1][i]=genome[mum][i];
						genome[genome.size()-2][i]=genome[dad][i];
					} 
					for(unsigned int i=0; i<genome[mum].size()-idx-1; i++){
						genome[genome.size()-1][idx+i] = mum_tail[i];
						genome[genome.size()-2][idx+i] = dad_tail[i];
					} 
					bool sent1 = check(genome[genome.size()-1]);
					bool sent2 = check(genome[genome.size()-2]);
					if(sent1){
						 cout << " Failure in crossover: mum will fix it " << endl << endl;
						 genome[genome.size()-1]=genome[mum];
					}
					if(sent2){
						 cout << " Failure in crossover: dad will fix it " << endl << endl;
						 genome[genome.size()-2]=genome[dad];
					}	 
					genome[genome.size()-1][genome[genome.size()-1].size()-1] = distance(genome[genome.size()-1]); 
					genome[genome.size()-2][genome[genome.size()-2].size()-1] = distance(genome[genome.size()-2]); 
				}
			}
		
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