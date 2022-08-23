#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

double OptimiseModel1 (run_params p, vector<int>& times, vector<double>& ratios, vector<double>& freqs, vector<treat>& med, gsl_rng *rgen);
void SetupInitialMu (run_params p, vector<double>& ratios, vector<double>& mu);
void CompareFitModel1 (run_params p, double& fit, double& best_fit, double& s, double& s_best, double& g, double& g_best, vector<double>& mu, vector<double>& mu_best);
void GetFuncParamsSimple (run_params p, const vector<double>& mu, const vector<treat>& med, vector<func_params_simple>& fp);
void FindOffsetsSimple (run_params p, const vector<double>& mu, const vector<treat>& med, vector<func_params_simple>& fp);
double CalcFreqPM (const double& u, const double& s, const double g, const double& o, const int& p, const int& t);
double CalcOffset (double q, double s, double g, double newu, double plus);
double CalculateFrequencySimple (run_params p, int t, const vector<treat>& med, const vector<func_params_simple>& fp);



