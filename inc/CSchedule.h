#ifndef __CSCH_H__
#define __CSCH_H__

#include <stdio.h>
#include <math.h>
class CState;


class CSchedule{
public:
	CSchedule();
	void SetTEnd(double tend){ tEnd = tend; }
	void SetdTmax(double dtmax){ dTmax = dtmax; }
	void SetReportTime(int _RepTimes, double tstep);
	void SetReportTime(int num_report_times, double * report_time);

	double GetTEnd(){ return tEnd; }
	double GetTStart(){ return tStart; }
	void AddStep(){ step++; }
	int GetStep(){ return step; }
	void SetTCurrent(double t){ tCur = t; }
	void SetTNext(double t){ tNext = t; }
	double GetDt(){ return dT; }
	void CalNewDt(CState *State, int converge);
	double GetTCurrent(){ return tCur; }
	double GetTNext(){ return tNext; }
	double GetCurrentReportTime(){ return currentReportTime; }
	int GetCurrentReportStep(){ return currentReportStep; }
	int get_num_report_times(){ return num_report_times_; }
	double get_report_time(int n){ return report_time_[n]; }
	void ChangedT(double dt){
		dTold = dt;
		dT = dt;
	}
	void ChangeOmege(double w){
		omega = w;
	}
	void ChangeYitaP(double np){
		yitaP = np;
	}
	void ChangeYitaS(double ns){
		yitaS = ns;
	}
private:
	int num_report_times_;
	double *report_time_;
	double currentReportTime;
	int currentReportStep;
	double dTold, dT, dTmax;
	double dateStart, dateEnd;
	double tStart, tEnd;
	double tCur, tNext;
	int step;

	double omega ;
	double yitaP ;
	double yitaS ;

};

#endif
