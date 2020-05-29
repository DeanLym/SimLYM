/*-----------------------------------------------------*/
/*  how to use the NOMAD library with a user function  */
/*-----------------------------------------------------*/
#include "sim_model.hpp"
#include "util_funs.hpp"


int main ( int argc , char ** argv ) {

	// display:
	int dim_ = 3600;
	if (argc < 2){
		cout << "Please specify input file." << std::endl;
	}
	ifstream ifs;
	ifs.open(argv[1]);
	vector<double> logk;
	double temp;
	for(int i=0;i<dim_;i++){
		ifs >> temp;
		logk.push_back(temp);
	}
	ifs.close();

	ifs.open(argv[2]);
	vector<int> well_location;
	int tmp;
	int num_well;
	int num_prod, num_inj;
	ifs >> num_prod;
	ifs >> num_inj;
	num_well = num_prod + num_inj;
	for(int i=0; i<2*num_well; i++){
		ifs >> tmp;
		well_location.push_back(tmp);
	}
	ifs.close();

	ifs.open(argv[3]);
	vector<double> bhp;
	for(int i=0; i<num_well; i++){
		ifs >> temp;
		bhp.push_back(temp);
	}
	ifs.close();

	vector<double> perm;
	perm.resize(dim_);
	vector<double> poro;
	poro.resize(dim_);
	GeneratePerm(dim_, &(logk[0]), &(perm[0]));

	SimCtrl* sim = GetSimulationModel(&perm[0], &well_location[0], &bhp[0], num_prod, num_inj);
	sim->display_level_ = 1;
	//cout << "Begin to run simulation." <<endl;
	sim->RunSim();
	//cout << "Finished running simulation." << endl;
	sim->OutputResult();
	return true;       // the evaluation succeeded
}
