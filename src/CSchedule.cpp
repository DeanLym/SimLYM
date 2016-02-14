

#include "CState.h"
#include <stdio.h>
#include <math.h>
#include"CSchedule.h"

CSchedule::CSchedule(){
	tStart = 0;
	dTold = 0.01;
	dT = 0.01;
	dTmax = 1e20;
	step = 0;
	omega = 0.5;
//	yitaP = 400;
	yitaP = 1000;
	yitaS = 0.40;
//	yitaS = 0.80;
}

void CSchedule::CalNewDt(CState *State, int converge){
	double min = 9999.9;
	double temp;
	const double *po = State->get_po();
	const double *po_n_ = State->get_po_n();
	const double *so = State->get_so();
	const double *so_n = State->get_so_n();
	if (converge == 0){
		dT = 0.5*dT;
	}
	else{
		dTold = dT;
		for (int i = 0; i < State->get_n_act_cell(); i++){
			temp = (1 + omega)*yitaP / (fabs(po[i] -po_n_[i]) + omega*yitaP);
			if (temp < min){
				min = temp;
			}
			temp = (1 + omega)*yitaS / (fabs(so[i] - so_n[i]) + omega*yitaS);
			if (temp < min){
				min = temp;
			}
		}
//		cout << min << endl;
		dT = dTold*min;
		if (dT > dTmax){
			dT = dTmax;
		}
		for (int i = 0; i<num_report_times_; i++){
			double flag_across_report_time;
			flag_across_report_time = (tNext - report_time_[i])*(tNext + dT - report_time_[i]);
			if (flag_across_report_time<0){
				dT = report_time_[i] - tNext;
				break;
			}
		}

	} // Only suitable for two phase problems
};

void CSchedule::SetReportTime(int num_report_times, double tstep){
	num_report_times_ = num_report_times;
	report_time_ = new double[num_report_times_];
	report_time_[0] = tstep;
	for (int i = 1; i < num_report_times_; i++){
		report_time_[i] = report_time_[i] + tstep;
	}
}

void CSchedule::SetReportTime(int num_report_times, double * report_time){
	num_report_times_ = num_report_times;
	report_time_ = new double[num_report_times_];
	for (int i = 0; i < num_report_times_; i++)
		report_time_[i] = report_time[i];
}
