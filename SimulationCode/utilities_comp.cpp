#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
using namespace std;
#include "simulation_run_infinite.h"
#include "utilities_comp.h"

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.u=0.1;
	p.s=0.1;
	p.g=400;
    p.jump=0;
    p.t_jump=30;
	p.mut_change=2;
	int x=1;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--u")==0) {
			x++;
			p.u=atof(argv[x]);
		} else if (p_switch.compare("--mut_change")==0) {
			x++;
			p.mut_change=atof(argv[x]);
        } else if (p_switch.compare("--t_jump")==0) {
            x++;
            p.mut_change=atof(argv[x]);
		} else if (p_switch.compare("--s")==0) {
			x++;
			p.s=atof(argv[x]);
        } else if (p_switch.compare("--jump")==0) {
            x++;
            p.jump=atoi(argv[x]);
		} else if (p_switch.compare("--g")==0) {
			x++;
			p.g=atoi(argv[x]);
		} else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void MakeInitialPopulation (run_params p, vector<double>& pop) {
    //The vector pop contains the fraction of genomes with pop[i] mutations
	for (int i=0;i<1001;i++) {
		pop.push_back(0);
	}
	pop[0]=1;
}

void SetupPoisson (double u, vector<double>& pois) {
    pois.clear();
	double tot=0;
    int limit=0;
    for (int i=0;i<100;i++) {
        tot=gsl_cdf_poisson_P(i,u);
        if (tot>0.999999) {
            limit=i;
            break;
        }
    }
    for (int i=0;i<=limit;i++) {
		double p=gsl_ran_poisson_pdf(i,u);
		pois.push_back(p);
	}
    tot=0;
    for (int i=limit+1;i<limit+10;i++) {
        double p=gsl_ran_poisson_pdf(i,u);
        tot=tot+(i*p);
    }
	pois[1]=pois[1]+tot;
}


//Next step - reimplement the mutation step using the multiple matrices
void MutationStep (run_params p, vector<double>& pois, vector<double>& pop) {
    //cout << "Mutation step\n";
    vector<double> oldpop=pop;
    for (int i=0;i<pop.size()-1;i++) { //e.g. pop.size()=50; from 0 to 49, so go to i=48.  At i=47 limit is 2
        int limit=min(pois.size(),pop.size()-i-1); //Maximum number of mutations to push
        //Amount coming out of population
        double total_out=0;
        for (int j=1;j<limit;j++) { //Calculate the amount of mutation going out of a point in the grid
            double am=pois[j]; //amount with j mutations
            total_out=total_out+am;
        }
        pop[i]=pop[i]-(total_out*oldpop[i]); //Remove things mutating out from pop[i][j]
        for (int j=1;j<limit;j++) {
            pop[i+j]=pop[i+j]+(pois[j]*oldpop[i]);
        }
    }
}

void FitnessStep (run_params p, vector<double>& pop) {
    double ftot=0;
    for (int i=0;i<pop.size();i++) {
        double fval=pow(p.s,i);
        ftot=ftot+(pop[i]*fval);
    }
    cout << "Total fitness " << ftot << "\n";
    for (int i=0;i<pop.size();i++) {
        double fval=pow(p.s,i);
        pop[i]=(pop[i]*fval)/ftot;
    }
}

double TotAllFreq (vector<double>& pop) {
    double tot=0;
    for (int i=0;i<pop.size();i++) {
        tot=tot+(i*pop[i]);
    }
    return tot;
}


void PrintPopulation(vector< vector<double> >& pop) {
	for (int i=0;i<pop[0].size();i++) {
		for (int j=0;j<pop.size();j++) {
			cout << pop[j][i] << " ";
		}
		cout << "\n";
	}
}
