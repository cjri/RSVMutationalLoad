#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>

void GetOptions (run_params& p, int argc, const char **argv);
void ReadFreqs (run_params p, vector<int>& times, vector<double>& freqs);
void ReadTreatment (run_params p, vector<treat>& med);

