#include "../libremia.h" 

bool comp(vector<double> pippo, vector<double> pluto){
		
			if( pippo[pippo.size()-1]<pluto[pluto.size()-1]) return true;
			else return false;
}


class Pop {
	
	private:
	
		vector<vector<double>> genome; //il genoma è un vettore di cromosomi (ognuno contenente una sequenza di città)
		vector<double> city; // sequenza contenente tutte le distanze delle città dalla prima
		double p_cross=0.5; // probabilità di crossing
		double p_mut=0.1; // probabilità di mutazione
		Random rand;
		bool dist;
		
	public:
		
		//Default constructor
		Pop(){}
		
		//constructor
		Pop(unsigned int n_chr, unsigned int n_genes){

			iniz(rand);
				
			vector<double> chromo(n_genes);
			city.push_back(0.); //fisso la prima città
			for(unsigned int i=0; i<n_genes; i++){
				chromo[i]=i;
				if(i!=0) city.push_back(rand.Rannyu(0.,2*M_PI)); //genero le altre città
			}
			
			sort(city.begin(),city.end()); //per fare i test, la sequenza corretta più breve tra le città è 1-2-3-4-5-6 ...
			
			for(unsigned int i=0; i<n_chr; i++){	//preparo i cromosomi con sequenze casuali
				random_shuffle(chromo.begin()+1, chromo.end());
				genome.push_back(chromo);
			}
			
		}
		
		//Destructor
		~Pop(){}
		
		//functions	
		
		int selection(){
			return int(genome.size()*pow(rand.Rannyu(),2.));
		}		
		
		void gensort(){
			sort(genome.begin(), genome.end(), comp); //uso il sort della standard library ma con una legge di scelta costruita ad hoc
		}
		
		
		//double eval1(const vector<double> chromo, const vector<double> city){
		double eval1(const vector<double> chromo){
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
			unsigned int appo = chromo.size()-1;
			index2 = chromo[appo];
			
			//if(city[index2]>M_PI) L += fabs(2*M_PI-city[index2]);
			//else L += fabs(city[index2]);
			
			Lx=cos(city[0])-cos(city[index2]);
			Ly=sin(city[0])-sin(city[index2]);
			L += sqrt((Lx*Lx)+(Ly*Ly));
			
			return L;
		}
		
		
		void mutations(){
			
			int chr = selection();
			
			if(rand.Rannyu()<p_mut){//scambio due città 
				int idx, idx2;
				idx = rand.Rannyu(1., genome[chr].size()-2);
				do{ idx2 = rand.Rannyu(1., genome[chr].size()-2); }while(idx==idx2);
				swap(genome[chr][idx], genome[chr][idx2]);
				genome[chr][genome[chr].size()-1] = eval1(genome[chr]);
				check(genome[chr]);
			}	
			/*
			chr = selection();
			
			if(rand.Rannyu()<p_mut){//sposto un blocco di 2 città di 4 posizioni
				int idx = rand.Rannyu(1., genome[chr].size()-8);
				int idx2 = idx + 2;
				vector<double> appo = {genome[chr][idx], genome[chr][idx+1]};
				for(unsigned int i=0; i<4; i++)	genome[chr][idx+i] = genome[chr][idx2+i];
				genome[chr][idx+4] = appo[0];
				genome[chr][idx+5] = appo[1];
				genome[chr][genome[chr].size()-1] = eval1(genome[chr]);
			}
			
			chr = selection();
			
			if(rand.Rannyu()<p_mut){//inverto l'ordine di 4 città
				int idx = rand.Rannyu(1., genome[chr].size()-2);
				vector<double> appo(4);
				for(unsigned int i=0; i<4; i++) appo[i] = genome[chr][idx+i];	
				for(unsigned int i=0; i<4; i++) genome[chr][idx+i] = appo[3-i];
				genome[chr][genome[chr].size()-1] = eval1(genome[chr]);
			
			}
			
			chr = selection();
			
			if(rand.Rannyu()<p_mut){//permuto l'ordine di 6 città
				int idx = rand.Rannyu(1., genome[chr].size()-7);
				random_shuffle(&genome[chr][idx], &genome[chr][idx+6]);
				genome[chr][genome[chr].size()-1] = eval1(genome[chr]);
			}
			*/
			
		}
		
		void crossover(){
			
			if(rand.Rannyu()<p_cross){
			
				int dad = selection();
				int mum = selection();
				
				int idx = rand.Rannyu(1, genome[mum].size()-2); //è un matriarcato
				
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
					check(genome[genome.size()-1]);
					check(genome[genome.size()-2]);
					genome[genome.size()-1][genome[genome.size()-1].size()-1] = eval1(genome[genome.size()-1]); 
					genome[genome.size()-2][genome[genome.size()-2].size()-1] = eval1(genome[genome.size()-2]); 
				}
			}
		
		}
		
		
		void check(vector<double> chr){
			double sum = 0;
			for(unsigned int i=0; i<chr.size()-1; i++) sum += chr[i];
			double realvalue = (city.size()-1)*city.size()/2;
			if(realvalue!=sum){ cerr << endl << "PROBLEM: chromosome is not correct" << endl;
			for(unsigned int i=0;  i<chr.size(); i++) cout << chr[i] << " ";
			cout << endl << endl;
			}
		}
		
		void setL1(){
			for(unsigned int i=0; i<genome.size(); i++){
				//genome[i].push_back(eval1(genome[i],city));
				genome[i].push_back(eval1(genome[i]));
			}
		}
		
		void set_probs(double p_c, double p_m){
			
			p_cross=p_c;
			p_mut=p_m;
		}
		
		void evolveL1(){
		
			int conta = 0;
			gensort();
		
			while(genome[0][genome[0].size()-1]>20){
				//mutations();
				crossover();
				gensort();
				
				conta++;
				
				//cout << conta << endl;
				if(conta%100==0){
				 cout << endl << " Ho appena generato le prime " << conta << " generazioni... "  << endl;
				 cout << endl << genome[0][genome[0].size()-1] << endl << endl;
				}
			}
		
		}
		
		void print_pop(){
		
			for(unsigned int j=0; j<genome.size(); j++){
			for(unsigned int i=0; i<genome[j].size(); i++) cout << genome[j][i] << " ";
			cout << endl << endl;
			}

		}
		
		void print_city(){
			for(unsigned int i=0; i<city.size(); i++) cout << i << " :  " << city[i] << endl;
			
			/*
			vector<double> chromo(city.size());
			for(unsigned int i=0; i<city.size(); i++) chromo[i]=i;
			cout << endl << endl << "L_min = " << eval1(chromo) << endl;
			
			/*
			Random rand;
			iniz(rand);
			for(unsigned int i=1; i<9*city.size(); i++) city[2]=rand.Rannyu(0.,2*M_PI);
			for(unsigned int i=1; i<city.size(); i++) city[i]=rand.Rannyu(0.,2*M_PI);
			
			sort(city.begin(),city.end());
			cout << endl << endl;
			for(unsigned int i=0; i<city.size(); i++) cout << i << " :  " << city[i] << endl;
			vector<double> chromo(city.size());
			for(unsigned int i=0; i<city.size(); i++) chromo[i]=i;
			double L=0.;
			unsigned int index1;
			unsigned int index2;
			
			for(unsigned int i=0; i<(chromo.size()-2); i++){
				index1 = chromo[i];
				index2 = chromo[i+1];
				L += fabs(city[index1]-city[index2]);
			}
			unsigned int appo = chromo.size()-1;
			index2 = chromo[appo];
			L += fabs(2.*M_PI-city[index2]);
			cout << endl << endl << "L_min = " << L << endl;
			*/
		}
};




int main(){

	Pop pop(100, 34);
	
	pop.setL1();
	pop.gensort();
	pop.print_pop();
	cout << endl << endl << endl << endl << endl;
	
	pop.evolveL1();
	//pop.gensort();
	pop.print_pop();
	/*for (unsigned int i=0; i<50; i++){
		pop.mutation();
		pop.print_pop();
	}*/
	
	//pop.print_city();

	return 0;
}
