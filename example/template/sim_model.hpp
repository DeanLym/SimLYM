/*
 * sim_model.hpp
 *
 *  Created on: Mar 5, 2016
 *      Author: yiminliu
 */

#ifndef SIM_MODEL_HPP_
#define SIM_MODEL_HPP_

#include "CSimCtrl.h"
#include <sstream>

using namespace::std;

SimCtrl* GetSimulationModel(double *kx, int *well_location, double* bhp, int num_prod, int num_inj){
	// Model 1 Simulator file
	SimCtrl* sim;
	sim = new SimCtrl;
	sim->display_level_ = 1;
	sim->state_report_ = 2;
	// input and initialize grid
	sim->grid_ = new CCartGrid(60,60,60);
	sim->grid_->InputDx(164.042);
	sim->grid_->InputDy(164.042);
	sim->grid_->InputDz(32.81);
	sim->grid_->InputTops(25574.15);
	sim->grid_->InputKx(kx);
	sim->grid_->InputKy(kx);
	sim->grid_->InputKz(kx);
	sim->grid_->InputPoro(0.2);
	sim->InitializeGrid();

	// initialize PVT
	sim->pvt_ = new CPVT;
	sim->pvt_->SetDensity(49.1, 64.79, 0.06);
	sim->pvt_->set_pvdo_table(42,"PVDO.DAT");
	sim->pvt_->set_pvtw_table("PVTW.DAT");

	double swi  = 0.10 , sor  = 0.30 ;
	double aw   = 1.5  , ao   = 3.6  ;
	double krw0 = 0.7  , kro0 = 1.0  ;
	sim->sat_   = new CSAT_COREY(swi, sor , aw ,ao ,krw0 ,kro0 );

	// initialize SCH
	sim->sch_ = new CSchedule;


	sim->sch_->SetTEnd(200.0);

	sim->sch_->SetTCurrent(0.0);
	sim->sch_->SetTNext(sim->sch_->GetDt() + sim->sch_->GetTCurrent());
	sim->sch_->SetdTmax(100.0);
	vector<double> report_time;
	int num_report_time = 1;

	num_report_time = 100;

	for(int i=0; i<num_report_time ;i++)
		report_time.push_back((i+1)*2);
	sim->sch_->SetReportTime(num_report_time, &report_time[0]);

	// initialize State
	sim->InitializeState();
	sim->SetInitPres(25590.55,4713.735);
//	sim->SetInitPres(po_file);
	// Set initial saturation equal to swi
	sim->SetInitSat(swi, 0.1);
//	sim->SetInitSw(sw_file);

	// initialize Well
	int ind_well_loc = 0;
	int ind_bhp = 0;
	for(int k=0; k<num_prod; k++){
		string well_name;
		stringstream ss;
		ss << "PROD-";
		ss << k + 1;
		ss >> well_name;

		CStandardWell *p_1 = new CSTDProdWell(well_name, sim->grid_->GetIndex(well_location[ind_well_loc],well_location[ind_well_loc+1],0));
		p_1->set_r(0.5);
		p_1->set_ctrl_mode(CStandardWell::CBHP);
		// cout << bhp[0] << endl;
		p_1->set_TL_BHP(bhp[ind_bhp]);
		p_1->CalWellIndex(sim->grid_);
		sim->std_well_.push_back(p_1);

		ind_well_loc += 2;
		ind_bhp += 1;
	}
	
	for(int k=0; k<num_inj; k++){
		string well_name;
		stringstream ss;
		ss << "INJ-";
		ss << k + 1;
		ss >> well_name;
		CStandardWell *inj_1 = new CSTDWInjWell(well_name, sim->grid_->GetIndex(well_location[ind_well_loc], well_location[ind_well_loc+1], 0));
		inj_1->set_r(0.5);
		inj_1->set_ctrl_mode(CStandardWell::CBHP);

		inj_1->set_TL_BHP(bhp[ind_bhp]);
		inj_1->CalWellIndex(sim->grid_);
		sim->std_well_.push_back(inj_1);

		ind_well_loc += 2;
		ind_bhp += 1;
	}

	//=====================================
	//	double wrat_target_inj = -9447.0, bhp_limit_inj = 14503.8;
	//	double bhp_inj = 7251.89;
	//	double bhp_inj = 5003.8;
	//	double bhp_inj = 5148.84;



	// initialize Solver
	sim->InitializeSolver();
	return sim;

}


#endif /* SIM_MODEL_HPP_ */
