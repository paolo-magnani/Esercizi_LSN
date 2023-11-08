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

void iniz(Random &rnd);

void iniz(Random &rnd, const unsigned int nrows);

double error(const vector<double> & AV,const vector<double> & AV2, unsigned int n);

void mediablocchi( const unsigned int M, const unsigned int N, const vector<double>& r, vector<double>& sum_prog, vector<double>& err_prog);

void stampamedia(string filename, const unsigned int& N_data, const vector<double>& media, const vector<double>& err);

vector<unsigned int> binsearch(const vector<double>& ran, const double& low, const double& high, const unsigned int& nbins);