#include<iostream>
#include<math.h>
#include<string>
#include<sstream>
#include <ctime>

using namespace std;
#include "model1.hpp"
#include "CSimCtrl.h"


int main(int argc, char ** argv){
	clock_t t1 = clock();
	char* perm_file = argv[1];
	string prod_bhp_str = argv[2];
	string inj_rate_str = argv[3];
	stringstream ss1;
	stringstream ss2;
	ss1 << prod_bhp_str;
	ss2 << inj_rate_str;
	double prod_bhp, inj_rate;
	ss1 >> prod_bhp;
	ss2 >> inj_rate;
	cout << prod_bhp << "  " << inj_rate << endl;
	SimCtrl* sim = model1(perm_file, prod_bhp, inj_rate);

	// Run Simulation
	clock_t t2 = clock();
	sim->RunSim();
	clock_t t3 = clock();
	sim->OutputResult();

	//
//	sim->hm_ = new CHistoryMatching;
//	sim->hm_->SetHMTarget("HIST.TXT");
//	double Sd = sim->hm_->GetDataMismatch(sim->std_well_);
//	cout << "Data mismatch is :" << Sd << endl;

	delete sim;
	clock_t t4 = clock();
	double t_init = double(t2 - t1) / CLOCKS_PER_SEC;
	double t_run = double(t3 - t2) / CLOCKS_PER_SEC;
	double t_total = double(t4 - t1) / CLOCKS_PER_SEC;
	double t_output = double(t4 - t3) / CLOCKS_PER_SEC;
	cout << endl << "Initilization Took " << t_init << " s" << endl;
	cout << endl << "Simulation Took " << t_run << " s" << endl;
	cout << endl << "Total Elpased Time:  " << t_total << " s" << endl;
	return 1;
}
