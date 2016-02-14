#include "CSimCtrl.h"
#include "CState.h"
#include "CPardisoSolver.h"
#include "CUnstructGrid.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace::std;

SimCtrl::SimCtrl()
:grid_type_(CARTESIAN),
 model_type_(SBLACKOIL),
 phase_type_(OW),
 implict_level_(FIM),
 linear_solver_(PARDISO),
 preconditioner_(NA),
 max_newton_iters_(10),
 min_err(1e-3),
 display_level_(1)
{

}

bool SimCtrl::InitializeGrid(){
	grid_->GenerateConnList();
	grid_->CalDepth();
	grid_->CalDeltaD();
	grid_->CalVolume();
	return true;
}

//bool SimCtrl::InitializeGrid(int ncell){
//	if (grid_type_ == UNSTRUCTURED){
//		grid_ = new CUnstructGrid(ncell);
//		return true;
//	}
//	else{
//		// Throwout an exception.
//		return false;
//	}
//}

//bool SimCtrl::InitializePVT(){
//	pvt_ = new CPVT;
//	pvt_->SetDensity(40.0, 62.238, 0.0);
//	return true;
//}

//bool SimCtrl::InitializeSAT(){
//	sat_ = new CSAT;
////	sat_->SetSWOF(40, "KR.INC");
//	return true;
//}

//bool SimCtrl::InitializeSCH(){
//	sch_ = new CSchedule;
//	sch_->SetTEnd(540.0);
//	sch_->SetTCurrent(0.0);
//	sch_->SetTNext(sch_->GetDt() + sch_->GetTCurrent());
//	sch_->SetReportTime(1, 540.0);
//	return true;
//}

bool SimCtrl::InitializeState(){
	state_ = new CState(phase_type_,grid_,pvt_);
	return true;
}

bool SimCtrl::SetInitPres(double datumDepth, double pDatum)
{
	state_->SetInitPres(datumDepth, pDatum, pvt_, grid_);
	return true;
}

bool SimCtrl::SetInitSat(double Sw, double Sg)
{
	state_->SetInitSat(Sw,Sg);
	return true;
}


bool SimCtrl::InitializeSolver(){
	if (linear_solver_ == PARDISO){
		solver_ = new CPardisoSolver(state_);
	}
	return true;
}



//bool SimCtrl::InitializeRun(){
//	set_grid_type(CARTESIAN);
//	set_phase_type(OW);
//	InitializePVT();
//	InitializeSAT();
//	InitializeSCH();
//	InitializeState();
//	InitializeWell();
//	InitializeSolver();
//	return true;
//}

bool SimCtrl::RunSim(){
	cout << endl << "Start of Simulation!" << endl << endl;
	double err;
	int n_iter; int temp;
	int converge;
	int honor_limit;
	const double *dx;
	double solver_time = 0.0;
	int count_newton = 0;
	state_->UpdateOldState(pvt_);
	while (sch_->GetTCurrent() < sch_->GetTEnd()){
		Display();
		//		cout << "Begin time step loop ..." << endl;
		converge = 1;
		n_iter = 0;
		err = min_err + 1;
		honor_limit = 1;
		while (err > min_err){
			//			cout << "Begin Newton loop ..." << endl;
			state_->CalParametersResidual(grid_, std_well_, pvt_, sat_, sch_);
			for (int i = 0; i < std_well_.size(); i++){
				//				cout << "Caluculating well residuals... " << endl;
				std_well_[i]->CalParameterResidual(state_);
			}
			err = state_->AssembleResidual(grid_, sch_->GetDt(), std_well_);
			//			cout << err << min_err << endl;
			//		    cout << "err:"<<err << endl;
			if (err < min_err){
				break;
			}
			state_->CalParametersJacobian(grid_, std_well_, pvt_, sat_, sch_);
			for (int i = 0; i < std_well_.size(); i++){
				std_well_[i]->CalParameterJacobian(state_);
			}
			state_->AssembleJacobian(grid_, sch_->GetDt(), std_well_);
			solver_->InputResidual(state_->get_resid());
			solver_->InputJacobian();
			//			cout << "Solving equation..." << endl;
			//			state_->OutputJacobian();
			dx = solver_->Solve();
			n_iter++;
			count_newton++;
			//			cout << "Newton Iteration " << n_iter << endl;
			if (n_iter > max_newton_iters_){
				converge = 0;
				break;
			}
			else{
				state_->UpdateNewState(dx);
			}
		}
		if (converge == 1){
			for (int i = 0; i < std_well_.size(); i++){
				if (!std_well_[i]->CheckLimits()){
					honor_limit = 0;
					//					cout << std_well_[i]->get_well_name() << " changed control mode" << endl;
				}
			}
			if (honor_limit){
				sch_->AddStep();
				double dt = sch_->GetDt();
				sch_->CalNewDt(state_, converge);
				sch_->SetTCurrent(sch_->GetTNext());
				sch_->SetTNext(sch_->GetTCurrent() + sch_->GetDt());
				state_->UpdateOldState(pvt_);
				for (int i = 0; i < std_well_.size(); i++){
					std_well_[i]->RecordResult(sch_->GetTCurrent());
				}
				tstep_.push_back(sch_->GetStep());
				num_newton_iter_.push_back(n_iter);
				cfl_.push_back(state_->CalCFLOW(grid_, dt));
				dt_.push_back(dt);


//				for (int i = 0; i < sch_->get_num_report_times(); i++){
//					if (sch_->GetTCurrent() == sch_->get_report_time(i)){
//						state_->State_Report(sch_, i + 1);
//						break;
//					}
//				}


				//Summary
				//				cout << "Time step:" << sch_->GetStep() << "  Current Time" << sch_->GetTCurrent() << endl << endl;
			}
			else{
				state_->ChangeBackState();
				converge = 0;
				sch_->CalNewDt(state_, converge);
				sch_->SetTNext(sch_->GetTCurrent() + sch_->GetDt());
			}
		}
		else{
			state_->ChangeBackState();
			sch_->CalNewDt(state_, converge);
			sch_->SetTNext(sch_->GetTCurrent() + sch_->GetDt());
		}
	}

	cout << endl << "End of Simulation!" << endl << endl;
	cout << "Total Number of Time Steps:  " << sch_->GetStep() << endl << endl;
	cout << "Total Number of Newton Iterations:  " << count_newton << endl << endl;

	return true;
}

bool SimCtrl::set_grid_type(GRIDTYPE grid_type){
	grid_type_ = grid_type;
	return true;
}

bool SimCtrl::set_phase_type(PHASETYPE phase_type){
	phase_type_ = phase_type;
	return true;
}

bool SimCtrl::OutputResult(){
	string fn;
	string sfx(".txt");
	ofstream out;
	vector<double> TIME, ORAT, WRAT, BHP,WBP;
	for (int i = 0; i < std_well_.size(); i++){
		fn = std_well_[i]->get_well_name();
		fn += sfx;
		out.open(fn.c_str());
		out << "TIME\t\t\t" << "ORAT\t\t\t" << "WRAT\t\t\t" << "BHP\t\t\t" << "Block Pressure" << endl;
		out << "DAT\t\t\t" << "STB/DAT\t\t\t" << "STB/DAT\t\t\t" << "PSI\t\t\t" << "PSI\t\t\t" <<endl;
		TIME = std_well_[i]->get_TIME();
		ORAT = std_well_[i]->get_ORAT();
		WRAT = std_well_[i]->get_WRAT();
		BHP = std_well_[i]->get_BHP();
		WBP = std_well_[i]->get_WBP();
		for (int k = 0; k < TIME.size(); k++){
			out << TIME[k] << "\t\t\t";
			out << ORAT[k] << "\t\t\t";
			out << WRAT[k] << "\t\t\t";
			out << BHP[k] << "\t\t\t";
			out << WBP[k] << "\t\t\t";
			out << endl;
		}
		out.close();
	}


	out.open("cfl.txt");
	out << "STEP\t\t\t" << "DT\t\t\t" << "CFL\t\t\t" << "NEWTON ITERATIONS\t\t\t" << endl;
	for (int i = 0; i < tstep_.size(); i++){
		out << tstep_[i] << "\t\t\t";
		out << dt_[i] << "\t\t\t";
		out << cfl_[i] << "\t\t\t";
		out << num_newton_iter_[i] << "\t\t\t";
		out << endl;
	}

	out.close();

	return true;
}

bool SimCtrl::Display(){
	if (display_level_ ==1)
		cout << "%% T = " <<sch_->GetTCurrent()
		<<"\t, DT = " << sch_->GetDt()<<endl;
	return true;
}
SimCtrl::~SimCtrl(){
	delete sat_;
	delete pvt_;
	delete grid_;
	delete sch_;
	delete state_;
	for (int i = 0; i < std_well_.size(); i++){
		delete std_well_[i];
	}
}
