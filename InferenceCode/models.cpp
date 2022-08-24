#include "basicmodel.h"
#include "models.h"
#include "utilities.h"
#include <string>

double OptimiseModel1 (run_params p, vector<int>& times, vector<double>& ratios, vector<double>& freqs, vector<treat>& med, gsl_rng *rgen) {
    //Optimise mu parameters and p.s
    
    //Set up initial mu
    vector<double> mu;
    SetupInitialMu (p,ratios,mu);
    
    //Set up for optimisation
    int first=1;
    double s=p.s;
    double s_best=s;
    double g=p.g;
    double g_best=g;
    vector<double> mu_best=mu;
    double dm=mu[0]/100;
    double ds=s/100;
    double dg=g/100;
	double best_fit=1e5;
	if (p.likelihood==1) {
		best_fit=-best_fit;
	}
	double fit=0;
    vector<func_params_simple> fp;

	for (int it=0;it<100000;it++) {
        if (first==0) {
            CompareFitModel1 (p,fit,best_fit,s,s_best,g,g_best,mu,mu_best);
        }
        first=0;

        //Change all parameters.
        s=s+(gsl_rng_uniform(rgen)*ds)-(ds/2);
        if (s>0.99) {
            s=0.99;
        }
       // g=g+(gsl_rng_uniform(rgen)*dg)-(dg/2);
        mu[0]=s*p.q_init;
        for (int m=1;m<mu.size();m++) {
            mu[m]=mu[m]+(gsl_rng_uniform(rgen)*dm)-(dm/2);
        }
     
        //Set up offsets
        p.s=s;
        p.g=g;
        GetFuncParamsSimple(p,mu,med,fp);
        FindOffsetsSimple(p,mu,med,fp);

       // cout << "Get freqs\n";
        fit=0;
        for (int i=0;i<times.size();i++) {
            double qn=CalculateFrequencySimple(p,times[i],med,fp);
            fit=fit+pow(freqs[i]-qn,2);
        }
        fit=sqrt(fit);
	}
    //Do something with the answer
    if (p.verb==1) {
        cout << "Best parameters: Mu ";
        for (int i=0;i<mu_best.size();i++) {
            cout << mu_best[i] << " ";
        }
        cout << "s " << s_best << " ";
        cout << "g " << g_best << " ";
        cout << " Fit " << best_fit << "\n";
        cout << "Frequencies\n";
    }
    p.s=s_best;
    mu=mu_best;
    //cout << "Size of mu is " << mu.size() << "\n";
    GetFuncParamsSimple(p,mu,med,fp);
    FindOffsetsSimple(p,mu,med,fp);
    if (p.verb==1) {
        for (int i=0;i<times.size();i++) {
            double qn=CalculateFrequencySimple(p,times[i],med,fp);
            cout << times[i] << " " << freqs[i] << " " << qn << "\n";
        }
    }
//    cout << "Best mu " << mu[0] << "\n";
//    cout << "Best s " << s_best << "\n";
//    cout << "Best g " << g_best << "\n";
    return s_best;
}

void SetupInitialMu (run_params p, vector<double>& ratios, vector<double>& mu) {
    mu.clear();
    cout << "Initial q " << p.q_init << "\n";
    p.mu=p.s*p.q_init; //Initial mu.
    mu.push_back(p.mu);
    for (int i=0;i<ratios.size();i++) {
        mu.push_back(ratios[i]*p.mu);
    }
}


void CompareFitModel1 (run_params p, double& fit, double& best_fit, double& s, double& s_best, double& g, double& g_best, vector<double>& mu, vector<double>& mu_best) {
    //RMSD calculation
    if (fit<best_fit) {
        best_fit=fit;
        mu_best=mu;
        s_best=s;
        g_best=g;
    } else {
        mu=mu_best;
        s=s_best;
        g=g_best;
    }
}

void GetFuncParamsSimple (run_params p, const vector<double>& mu, const vector<treat>& med, vector<func_params_simple>& fp) {
    //Single selection coefficient
    fp.clear();
    for (int i=1;i<med.size();i++) { //For each block of time: Potentially different mutation rate
        func_params_simple f;
        f.index=i;
        f.s=p.s;
        f.g=p.g/24;
        f.mu=mu[med[i].type];
        if (f.mu>mu[i-1]){
            f.plus=-1;
        } else {
            f.plus=1;
        }
        fp.push_back(f);
    }
}

void FindOffsetsSimple (run_params p, const vector<double>& mu, const vector<treat>& med, vector<func_params_simple>& fp) {
    for (int i=0;i<fp.size();i++) {
        fp[i].o=0;
        if (i==0) {
            //Based upon converged frequencies
            double qi=mu[0]/fp[i].s;
            /*if (p.verb==1) {
                cout << qi << " " << fp[i].s << " " << fp[i].mu << " " << fp[i].plus << "\n";
            }*/
            if (qi*fp[i].s/fp[i].mu>1) {
                fp[i].plus=1;
            } else {
                fp[i].plus=-1;
            }
            double off=CalcOffset(qi,fp[i].s,fp[i].g,fp[i].mu,fp[i].plus);
            fp[i].o=off-med[1].date;
        } else {
            //Based upon last frequency of previous treatment regime
            double qi=CalcFreqPM(fp[i-1].mu,fp[i].s,fp[i].g,fp[i-1].o,fp[i-1].plus,med[i+1].date);
            /*if (p.verb==1) {
                cout << qi << " " << fp[i].s << " " << fp[i].mu << " " << fp[i].plus << "\n";
            }*/
            if (qi*fp[i].s/fp[i].mu>1) {
                fp[i].plus=1;
            } else {
                fp[i].plus=-1;
            }

            double off=CalcOffset(qi,fp[i].s,fp[i].g,fp[i].mu,fp[i].plus);
            fp[i].o=off-med[i+1].date;
        }
        /*if (p.verb==1) {
            cout << "Offset " << fp[i].o << "\n";
        }*/
    }
}

double CalcFreqPM (const double& u, const double& s, const double g, const double& o, const int& p, const int& t) {
    double q1=u/s;
    double q2=(1+(p*exp(-(s/g)*(t+o))));
    double q=q1*q2;
   // cout << q1 << " " << q2 << " " << t+o << " " << q << "\n";
    return q;
}


double CalcOffset (double q, double s, double g, double newu, double plus) {
    double o=q*s/newu;
    //cout << "O " << q << " " << s << " " << newu << " " << o << " " << plus;
    if (abs(o-1)<1e-6) {
        o=0;
    } else if (plus==-1) {
        o=log(1-o);
    } else {
        o=log(o-1);
    }
    o=-o/(s/g);
    //cout << " " << o << "\n";
    return o;
}

double CalculateFrequencySimple (run_params p, int t, const vector<treat>& med, const vector<func_params_simple>& fp) {
    double qn=0;
    if (t<med[1].date) {
        qn=p.q_init;
    } else {
        for (int i=1;i<med.size();i++) {
//            cout << "Now here " << fp[i-1].mu << " s " << fp[i-1].s << " o " << fp[i-1].o << " p " << fp[i-1].plus << "\n";
            int max=10000;
            if (i<med.size()-1) {
                max=med[i+1].date;
            }
            if (t>=med[i].date&&t<max) {
                qn=CalcFreqPM(fp[i-1].mu,fp[i-1].s,fp[i-1].g,fp[i-1].o,fp[i-1].plus,t);
                break;
            }
        }
    }
    //cout << "Qn is " << qn << "\n";
    return qn;
}


//Add in a module here to fit a likelihood to the inference/data combination











