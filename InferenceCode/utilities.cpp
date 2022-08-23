#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>

#include "basicmodel.h"
#include "utilities.h"
#include "io.h"

void GetWindowEnds (const vector<treat>& med, const vector<int>& times, vector<int>& window_ends) {
    //Treatment windows
    for (int i=1;i<med.size();i++) {
        for (int j=0;j<times.size();j++) {
            if (times[j]>med[i].date) {
                window_ends.push_back(j-1);
                break;
            }
        }
    }
    window_ends.push_back(times.size()-1);
}

double CalculateQInit (const vector<double>& freqs, const vector<int> window_ends, vector<double>& initial) {
	double q=0;
	for (int i=0;i<=window_ends[0];i++) {
		q=q+freqs[i];
        initial.push_back(freqs[i]);
	}
	q=q/(window_ends[0]+1);
	return(q);
}

void FindSigma (run_params p, double& sigma, const vector<double>& initial) {
    //Calculate standard deviation from the initial data
    for (int i=0;i<initial.size();i++) {
        //cout << p.q_init << " " << initial[i] << "\n";
        sigma=sigma+pow((p.q_init-initial[i]),2);
    }
    sigma=sigma/(initial.size()+0.);
    sigma=sqrt(sigma);
}

void FindRatios (run_params p, int n_types, const vector<double>& freqs, const vector<int> window_ends, vector<treat> med, vector<double>& ratios) {
    //Find mean frequency towards the end of a window
    //Use this to set the initial ratio between mutation rates
    for (int w=1;w<window_ends.size();w++) {
        int start=window_ends[w-1]+1;
        int end=window_ends[w];
        int half=start+((end-start)/2);
        double q=0;
        double count=0;
        for (int i=half;i<=end;i++) {
            q=q+freqs[i];
            count++;
        }
        if (count==0) {
            ratios.push_back(1);
        } else {
            q=q/count;
            ratios.push_back(q/p.q_init);
        }
    }
	
	//Additional code: Modify ratios according to the number of treatments
	
	//Set up vector for new ratios
	vector<double> new_ratios;
	vector<double> nrcount;
	for (int i=0;i<n_types;i++) {
		new_ratios.push_back(0);
		nrcount.push_back(0);
	}
	
	//Mean of intial ratios assigned for each type of treatment
	for (int i=0;i<ratios.size();i++) {
		new_ratios[med[i+1].type-1]=new_ratios[med[i+1].type-1]+ratios[i];
		//cout << i << " " << med[i+1].type << " " << new_ratios[med[i+1].type-1] << "\n";
		nrcount[med[i+1].type-1]++;
	}
	for (int i=0;i<new_ratios.size();i++) {
		//cout << "i " << i << " " << new_ratios[i] << " " << nrcount[i] << "\n";
		new_ratios[i]=new_ratios[i]/nrcount[i];
	}

	ratios=new_ratios;
	
}

void FindMuSize(const vector<treat>& med, int& mu_size) {
    for (int i=0;i<med.size();i++) {
        if (med[i].type>mu_size) {
            mu_size=med[i].type;
        }
    }
}


double CalculateLikelihood (const vector<int>& times, const vector<double>& freqs, const vector<double>& func) {
	double logL=0;
	for (int i=0;i<freqs.size();i++) {
		//cout << times[i] << " " << freqs[i] << " " << func[times[i]] << "\n";
		//logL=logL+gsl_ran_gaussian_pdf(freqs[i]-func[times[i]],s);
		logL=logL-pow(freqs[i]-func[times[i]],2);
	}
	return logL;
}
