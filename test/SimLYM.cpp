#include<iostream>
#include<math.h>
#include<string>
#include <ctime>

using namespace std;
#include "CSimCtrl.h"


int main(){
	clock_t t1 = clock();
	SimCtrl* sim;
	sim = new SimCtrl;
	// input and initialize grid
	sim->grid_ = new CCartGrid(60,60,1);
	sim->grid_->InputDx(82.02);
	sim->grid_->InputDy(82.02);
	sim->grid_->InputDz(32.81);
	sim->grid_->InputTops(11466.54);
	sim->grid_->InputKx("PERM.DAT");
	sim->grid_->InputKy("PERM.DAT");
	sim->grid_->InputKz("PERM.DAT");
	sim->grid_->InputPoro(0.2);
	sim->InitializeGrid();
	sim->grid_->OutputConnList();

	// initialize PVT
	sim->pvt_ = new CPVT;
	sim->pvt_->SetDensity(49.1, 64.79, 0.06);
	sim->pvt_->set_pvdo_table(42,"PVDO.DAT");
	sim->pvt_->set_pvtw_table("PVTW.DAT");

	// initialize SAT
	sim->sat_ = new CSAT;
	sim->sat_->SetSWOF(9,"SCAL.DAT");

	// initialize SCH
	sim->sch_ = new CSchedule;
	sim->sch_->SetTEnd(540.0);
	sim->sch_->SetTCurrent(0.0);
	sim->sch_->SetTNext(sim->sch_->GetDt() + sim->sch_->GetTCurrent());
	sim->sch_->SetdTmax(30.0);
	vector<double> report_time;
	report_time.push_back(1); report_time.push_back(2);report_time.push_back(5);
	sim->sch_->SetReportTime(3, &report_time[0]);

	// initialize State
	sim->InitializeState();
	sim->SetInitPres(11482.94,5076.33);
	sim->SetInitSat(0.1, 0.0);


	// initialize Well
	CStandardWell *p_1 = new CSTDProdWell("PROD-1", sim->grid_->GetIndex(4,19,0));
	p_1->set_r(1);
	p_1->set_ctrl_mode(CStandardWell::CBHP);
	p_1->set_TL_BHP(2900.76);
	p_1->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_1);
	CStandardWell *p_2 = new CSTDProdWell("PROD-2", sim->grid_->GetIndex(9,9,0));
	p_2->set_r(1);
	p_2->set_ctrl_mode(CStandardWell::CBHP);
	p_2->set_TL_BHP(2900.76);
	p_2->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_2);
	CStandardWell *p_3 = new CSTDProdWell("PROD-3", sim->grid_->GetIndex(54,54,0));
	p_3->set_r(1);
	p_3->set_ctrl_mode(CStandardWell::CBHP);
	p_3->set_TL_BHP(2900.76);
	p_3->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_3);
	//=====================================
	CStandardWell *inj_1 = new CSTDWInjWell("INJ-1", sim->grid_->GetIndex(29,24, 0));
	inj_1->set_r(1);
	inj_1->set_ctrl_mode(CStandardWell::CBHP);
	inj_1->set_TL_BHP(7251.9);
	inj_1->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(inj_1);

	CStandardWell *inj_2 = new CSTDWInjWell("INJ-2", sim->grid_->GetIndex(29,44, 0));
	inj_2->set_r(1);
	inj_2->set_ctrl_mode(CStandardWell::CBHP);
	inj_2->set_TL_BHP(8702.28);
	inj_2->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(inj_2);

	// initialize Solver
	sim->InitializeSolver();

	// Run Simulation
	clock_t t2 = clock();
	sim->RunSim();
	clock_t t3 = clock();
	sim->OutputResult();
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
