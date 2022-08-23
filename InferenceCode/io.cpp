#include "basicmodel.h"
#include "io.h"
#include <string>


void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.mu=1e-5;
	p.g=12;
	p.s=0.1;
    p.freqs_file="Freqs.dat";
    p.treat_file="Treatment.dat";
	p.likelihood=0;
	int x=1;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		} else if (p_switch.compare("--g")==0) {
			x++;
			p.g=atof(argv[x]);
		} else if (p_switch.compare("--s")==0) {
			x++;
			p.s=atof(argv[x]);
        } else if (p_switch.compare("--freqs")==0) {
            x++;
            p.freqs_file=argv[x];
        } else if (p_switch.compare("--treat")==0) {
            x++;
            p.treat_file=argv[x];
        } else if (p_switch.compare("--mu2")==0) {
            x++;
            p.mu2=atof(argv[x]);
       } else if (p_switch.compare("--check")==0) {
            x++;
            p.check=atoi(argv[x]);
       } else if (p_switch.compare("--like")==0) {
            x++;
            p.likelihood=atoi(argv[x]);
       } else if (p_switch.compare("--sigma")==0) {
            x++;
            p.sigma=atof(argv[x]);
		} else {
			cout << "Incorrect usage " << argv[x] << "\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void ReadFreqs (run_params p, vector<int>& times, vector<double>& freqs) {
	ifstream freqs_file;
	freqs_file.open(p.freqs_file);
	int n;
	double x;
	for (int i=0;i<100000;i++) {
		if (!(freqs_file >> n)) break;
		if (!(freqs_file >> x)) break;
		times.push_back(n);
		freqs.push_back(x);
	}
}

void ReadTreatment (run_params p, vector<treat>& med) {
	ifstream freqs_file;
	freqs_file.open(p.treat_file);
	int n;
	for (int i=0;i<100000;i++) {
		treat t;
		if (!(freqs_file >> n)) break;
		t.date=n;
		if (!(freqs_file >> n)) break;
		t.type=n;
		med.push_back(t);
	}
}
