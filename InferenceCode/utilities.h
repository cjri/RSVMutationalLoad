#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

void GetWindowEnds (const vector<treat>& med, const vector<int>& times, vector<int>& window_ends);
double CalculateQInit (const vector<double>& freqs, const vector<int> window_ends, vector<double>& initial);
void FindSigma (run_params p, double& sigma, const vector<double>& initial);
void FindRatios (run_params p, int n_types, const vector<double>& freqs, const vector<int> window_ends, vector<treat> med, vector<double>& ratios);
void FindMuSize(const vector<treat>& med, int& mu_size);
double CalculateLikelihood (const vector<int>& times, const vector<double>& freqs, const vector<double>& func);
