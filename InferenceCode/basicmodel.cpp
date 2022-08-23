#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>


#include "basicmodel.h"
#include "models.h"
#include "utilities.h"
#include "io.h"

int main(int argc, const char **argv) {

	int seed=(int) time(NULL);
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);
	
	run_params p;
	GetOptions (p,argc,argv);
	
	//Read in frequency data
	vector<int> times;
	vector<double> freqs;
	ReadFreqs(p,times,freqs);
    
	//Read in times of treatment
	vector<treat> med;
	ReadTreatment(p,med);
	int n_types=0;
	for (int i=0;i<med.size();i++) {
		if (med[i].type>n_types) {
			n_types=med[i].type;
		}
	}

	//Build list of samples within each treatment window
	vector<int> window_ends;
    GetWindowEnds(med,times,window_ends);
	
	//Calculate initial mean frequency - won't change at all
    //N.B. Require at least one sample with no treatment to validate
    vector<double> initial;
    double q=CalculateQInit(freqs,window_ends,initial);
    p.q_init=q;
    double sigma=0;
    FindSigma(p,sigma,initial); //Useful if likelihood needed
   
    //Find ratios to guess initial mu values.  These initialise the mutation rates used in the optimisation.
    vector<double> ratios;
	FindRatios (p,n_types,freqs,window_ends,med,ratios);
	
    //Note: Code may lose generality for other datasets e.g. if there are fewer data points.
    
    /*cout << "Ratios\n";
    for (int i=0;i<ratios.size();i++) {
        cout << ratios[i] << "\n";
    }*/
	
    int mu_size=0;
    FindMuSize(med,mu_size);
  
    //Fit a deterministic model
    
    
    //Use a pre-determined fitness value s to learn initial u.
    p.verb=1;
    OptimiseModel1 (p,times,ratios,freqs,med,rgen);

    
    return 0;
}
