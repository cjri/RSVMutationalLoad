using namespace std;
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>


struct run_params {
	double mu; //Initial mutatation rate
	double g; //Generation time
	double s; //Cost of mutation
	double q_init; //Initial frequency
    string freqs_file; //Name of file with frequency data
    string treat_file; //Name of file with treatment data
    int check; //Function to check an output
    double mu2; //For input of next mu value if required
    int likelihood; //Fit to data using normal distribution
    double sigma; //Standard deviation for likelihood model
    int verb; //Verbose output
};

struct func_params_simple {
    int index;
    double mu;
    double o;
    double s;
    double g;
    int plus;
};

	
struct func_params {
    int index;
    double mu;
    vector<double> o;
    vector<double> s;
    vector<int> plus;
};

struct treat {
	int date;
	int type;
};
