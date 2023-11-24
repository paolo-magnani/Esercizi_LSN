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

void iniz(Random &rnd); //inizializza il generatore di numeri casuali con i due numeri primi del file Primes della prima riga

void iniz(Random &rnd, const unsigned int nrows); //inizializza il generatore di numeri casuali con due numeri primi del file Primes posizionati nella riga nrows (da 1 a 384)

double error(const vector<double> & averages,const vector<double> & squaredaverages, unsigned int ndata); //calcola la deviazione standard della media su un blocco

void mediablocchi(const unsigned int& M, const unsigned int& N, const vector<double>& data, vector<double>& mean, vector<double>& err); //calcola la media a blocchi su M dati, suddividendoli in N blocchi e immagazina le medie e gli errori in mean e err

void stampamedia(string filename, const unsigned int& M, const vector<double>& media, const vector<double>& err); //stampa su file l'andamento progressivo della media a blocchi sugli M dati

vector<unsigned int> binsearch(const vector<double>& data, const double& min_value, const double& max_value, const unsigned int& nbins); // algoritmo di ricerca binaria su dati compresi tra un minimo e un massimo e incasellati in nbins