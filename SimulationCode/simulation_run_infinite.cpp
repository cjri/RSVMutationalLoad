//Simulate the effect of a change in mutation rate on allele frequency

#include <iostream>
#include <vector>
#include <list>
#include <string>
using namespace std;

#include "simulation_run_infinite.h"
#include "utilities_comp.h"

/*#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>*/

vector<double> p_select;

int main(int argc, const char **argv){

	run_params p;
	GetOptions(p,argc,argv);

    //Set up random number generator
	int seed=time(NULL);
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);

    double u=p.u; //Mutation rate for the whole genome per generation i.e. mu * L
    vector<double> pois;
    //Have a vector that describes the proportion of individuals getting 1, 2, 3 mutations etc.
    
    //Find maximum Poisson size
    SetupPoisson(u,pois);

    //Setup initial population
    vector<double> pop;
    MakeInitialPopulation(p,pop);
 
    /*cout << "Poisson ";
    for (int i=0;i<pois.size();i++) {
        cout << pois[i] << " ";
    }
    cout << "\n";*/
    
    double tot=0; //Total mutational load: Sum of allele frequencies
   // double tot_prev=0;  //Can be used to detect convergence - see commented out section below

    vector<double> freqs;
    vector<double> pop_sample=pop;

	for (int gen=0;gen<p.g;gen++) {
            //Mutation
            MutationStep (p,pois,pop);
            //Calculate total allele frequency
            tot=TotAllFreq(pop);
            cout << "Sum allele frequency " << tot << " ";
            
            freqs.push_back(tot);
            /*if (abs(tot_prev-tot)<(tot*1e-7)) {
                FitnessStep(p,pop);
                break;
            } else {
                tot_prev=tot;
            }*/
        //Fitness step
        FitnessStep(p,pop);
	}

    //Jump step: Following equilibrium, model 200 generations of evolution
    if (p.jump==1) {
        for (int gen=0;gen<p.g;gen++) {
            if (gen==p.t_jump) { //Code to change mutation rate
                u=u*p.mut_change;
                pois.clear();
                SetupPoisson(u,pois);
                /*cout << "Poisson ";
                for (int i=0;i<pois.size();i++) {
                    cout << pois[i] << " ";
                }
                cout << "\n";*/
            }
            MutationStep (p,pois,pop);
            //Calculate total allele frequency
            tot=TotAllFreq(pop);
            cout << "Sum allele frequency " << tot << " ";
            freqs.push_back(tot);
            //Fitness step
            FitnessStep(p,pop);
        }
    }
    
	return 0;

}
	
	
