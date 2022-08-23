
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>

struct run_params {
    int g; //Number of generations to simulate
    double u; //Mutation rate per genome i.e. /mu * L
    double s; //Constant fitness cost
    int jump; //Change mutation rate
    int t_jump; //Time of change in mutation rate
	double mut_change; //Multiplier for mutation rate
};

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>

void GetOptions (run_params& p, int argc, const char **argv);
void MakeInitialPopulation (run_params p, vector<double>& pop);
void SetupPoisson (double u, vector<double>& pois);
void MutationStep (run_params p, vector<double>& pois, vector<double>& pop);
void FitnessStep (run_params p, vector<double>& pop);
double TotAllFreq (vector<double>& pop);
void PrintPopulation(vector< vector<double> >& pop);




