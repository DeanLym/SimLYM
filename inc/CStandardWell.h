#ifndef __CSTANDARDWELL_H__
#define __CSTANDARDWELL_H__

#include <string>
#include <vector>
using namespace std;

class CGrid;
class CState;

class CStandardWell{
public:
	enum CTRLMODE{CBHP,CORAT,CWRAT,CLRAT,CGRAT};
	enum STDWELLTYPE{STDPROD,STDWINJ};
	static const double PI = 3.141592653;
	static const double alpha = 0.001127;
public:
	CStandardWell(string well_name,int block_index);
	~CStandardWell();
public:
	string get_well_name() const;
	bool set_r(double r);
	bool CalWellIndex(CGrid* grid);
	bool set_well_index(double WI);
	int get_n();
public:
	bool set_ctrl_mode(CTRLMODE ctrl_mode);
	bool set_TL_ORAT(double TL_ORAT);
	bool set_TL_BHP(double TL_BHP);
	bool set_TL_WRAT(double TL_WRAT);
	bool set_TL_GRAT(double TL_GRAT);
	bool set_TL_LRAT(double TL_LRAT);
public:
	double get_qo();
	double get_qw();
	double get_qg();
	double get_dqo_dp();
	double get_dqw_dp();
	double get_dqg_dp();
	double get_dqo_dsw();
	double get_dqw_dsw();
	double get_dqg_dsw();
	double get_dqo_dsg();
	double get_dqw_dsg();
	double get_dqg_dsg();
public:
	virtual bool CalParameterResidual(CState *State) = 0;
	virtual bool CalParameterJacobian(CState *State) = 0;
	virtual bool CalBHP(CState *state) = 0;
	virtual bool CalORAT(CState *state) = 0;
	virtual bool CalWRAT(CState *state) = 0;
	virtual bool CalGRAT(CState *state) = 0;
	virtual bool CalLRAT(CState *state) = 0;
	virtual bool CalDqoDp(CState *State) = 0;
	virtual bool CalDqwDp(CState *State) = 0;
	virtual bool CalDqgDp(CState *State) = 0;
	virtual bool CalDqoDsw(CState *State) = 0;
	virtual bool CalDqwDsw(CState *State) = 0;
public:
	virtual bool CheckLimits() = 0;
public:
	bool RecordResult(double t);
public:
	vector<double> get_TIME() const;
	vector<double> get_ORAT() const;
	vector<double> get_WRAT() const;
	vector<double> get_BHP() const;
	vector<double> get_LRAT() const;
	vector<double> get_GRAT() const;
	vector<double> get_WBP() const;
	vector<double> get_WWCT() const;
protected:
	string well_name_;
	STDWELLTYPE well_type_;
	int n_;
	//=================================
	double WI_;
	double r_, r0_;

	double po_, pg_, pw_;
	double mo_, mg_, mw_;
	double bo_, bg_, bw_;
	double rso_;
	double dmo_dp_, dmg_dp_, dmw_dp_;
	double dmo_dsw_, dmw_dsw_, dmg_dsw_;
	double dmo_dsg_, dmg_dsg_;

	double ORAT_, GRAT_, WRAT_, BHP_, LRAT_;
	double ORAT_RC_, LRAT_RC_, WRAT_RC_, GRAT_RC_;
	double dqo_dp_, dqo_dsw_, dqo_dsg_;
	double dqw_dp_, dqw_dsw_, dqw_dsg_;
	double dqg_dp_, dqg_dsw_, dqg_dsg_;

	//=================================
	CTRLMODE ctrl_mode_;
	double TL_ORAT_, TL_WRAT_, TL_GRAT_, TL_BHP_, TL_LRAT_;


	//Result vector
	vector<double> TIME_;
	vector<double> ORAT_REC_;
	vector<double> ORAT_RC_REC_;
	vector<double> GRAT_REC_;
	vector<double> GRAT_RC_REC_;
	vector<double> WRAT_REC_;
	vector<double> WRAT_RC_REC_;
	vector<double> BHP_REC_;
	vector<double> LRAT_REC_;
	vector<double> LRAT_RC_REC_;
	vector<double> WWCT_REC_;
	vector<double> WBP_REC_; // Well Block Pressure


};

#endif
