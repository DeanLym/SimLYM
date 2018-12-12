#ifndef __CState_H__
#define __CState_H__

#include"CPVT.h"
#include"CSAT.h"
#include"CGrid.h"
#include"CMatrixBlock.h"
#include"CConn.h"
#include"CSimCtrl.h"
#include <iostream>
#include <vector>
#include "CMatrixBlock.h"

#include "H5Cpp.h"

using namespace H5;

using namespace std;

class CJacob;
class CSchedule;
class MatrixBlock;
class CStandardWell;

class CState{
public:
	static const double g_ = 32.2;
	static const double gc = 32.2;
	static const double beta = 0.0069444444;
	static const double M = 5.615;
public:
	CState(SimCtrl::PHASETYPE phase, CGrid *Grid, CPVT* PVT);
	~CState();
	int get_nconns(){ return nconns_; }
	int get_num_phase(){ return num_phase_; }
	int GetHasOil(){ return hasOil; }
	int GetHasWater(){ return hasWater; }
	int GetHasGas(){ return hasGas; }
	int get_n_act_cell(){ return n_act_cell_; }

	bool SetInitPres(double datumDepth, double pDatum, CPVT *PVT, CGrid *Grid);
	bool SetInitPres(double po);
	bool SetInitPres(double *po);
	bool SetInitPres(char *po_file);

	bool SetInitSat(double Sw, double Sg);
	bool SetInitSat(double *Sw, double *Sg);
	bool SetInitSat(char *Sw_file, char *Sg_file);
	bool SetInitSw(char *Sw_file);
	bool SetInitSg(char *Sg_file);
	
	const double* GetBo(){ return bo_; }
	const double* GetBg(){ return bg_; }
	const double* GetBw(){ return bw_; }

    const double* GetKro(){ return kro_; }
    const double* GetKrg(){ return krg_; }
    const double* GetKrw(){ return krw_; }

    const double* GetVisO(){ return vo_; }
	const double* GetVisG(){ return vg_; }
	const double* GetVisW(){ return vw_; }

    bool CalDPhiOil(CGrid *Grid);
    bool CalDPhiGas(CGrid *Grid);
    bool CalDPhiWater(CGrid *Grid);
    bool CalGammaOil();
    bool CalGammaGas();
    bool CalGammaWater();

    bool CalMobiWater();
    bool CalAccuOil(CGrid *Grid, double dt);
    const double* GetAccuOil(){ return ao_; }
    bool CalAccuWater(CGrid *Grid, double dt);
    const double* GetAccuWater(){ return aw_; }
	bool CalAccuGas(CGrid* Grid, double dt);
    const double* GetAccuGas(){ return ag_; }

	const double* get_so(){ return so_; }
	const double* get_so_n(){ return so_n_; }
	const double* get_sw(){ return sw_; }
	const double* get_sw_n(){ return sg_n_; }
	const double* get_sg(){ return sg_; }
	const double* get_sg_n(){ return sg_n_; }
	const double* get_po(){ return po_; }
	const double* get_po_n(){ return po_n_; }

	double get_bo(int n){ return bo_[n]; }
	double get_bw(int n){ return bw_[n]; }
	double get_bg(int n){ return bg_[n]; }
	double get_rso(int n){ return rso_[n]; }

	double get_mo(int n){ return mo_[n]; }
	double get_mw(int n){ return mw_[n]; }
	double get_mg(int n){ return mg_[n]; }

	double get_dmo_dpo(int n){ return dmo_dpo_[n]; }
	double get_dmw_dpw(int n){ return dmw_dpw_[n]; }
	double get_dmg_dpg(int n){ return dmg_dpg_[n]; }
	double get_dmo_dsw(int n){ return dmo_dsw_[n]; }
	double get_dmo_dsg(int n){ return dmo_dsg_[n]; }
	double get_dmw_dsw(int n){ return dmw_dsw_[n]; }
	double get_dmg_dsg(int n){ return dmg_dsg_[n]; }



	const double* get_pg(){ return pg_; }
	const double* get_pg_n(){ return pg_n_; }
	const double* get_pw(){ return pw_; }
	const double* get_pw_n(){ return pw_n_; }

    bool CalMobiOil();
    bool CalMobiGas();

	void UpdateNewState(const double* dx);
    bool UpdateOldState(CPVT *PVT);
	void ChangeBackState();

    void Init_State_Report(int state_report);
	void State_Report(CSchedule *SCH,int n);
	void State_Debug();
	void Var_Debug(double *var, int n, char* fn);
	void Residual_Debug(char* fn);
	void CalParametersResidual(CGrid *Grid, vector<CStandardWell*> &std_well, CPVT *PVT, CSAT *SAT, CSchedule *SCH);
	void CalParametersJacobian(CGrid *Grid, vector<CStandardWell*> &std_well, CPVT *PVT, CSAT *SAT, CSchedule *SCH);
    void OutPutFPR(CSchedule &SCH);

    bool CalFo();
    bool CalFw();
    bool CalFg();

	bool CalDMoDPo();
	bool CalDMwDPw();

	bool CalDMoDSw();
	bool CalDMwDSw();

	bool CalDGammaODPo();
	bool CalDGammaWDPw();

	bool CalDPoroDPo(CGrid *Grid);

	double CalResidualNormInf(CGrid *Grid, double dt);

	double AssembleResidual(CGrid* Grid, double dt, vector<CStandardWell*> &std_well);
	bool AssembleJacobian(CGrid* Grid, double dt, vector<CStandardWell*> &std_well);

	bool AssembleResidualOG(vector<CStandardWell*> &std_well);
	bool AssembleResidualOW(vector<CStandardWell*> &std_well);
	bool AssembleResidualOWG(vector<CStandardWell*> &std_well);

	bool AssembleJacobianOG(CGrid* Grid, double dt, vector<CStandardWell*> &std_well);
	bool AssembleJacobianOW(CGrid* Grid, double dt, vector<CStandardWell*> &std_well);
	bool AssembleJacobianOWG(CGrid* Grid, double dt, vector<CStandardWell*> &std_well);

	double CalCFLOW(CGrid* Grid, double dt);

	const int* get_l(){ return l_; }
	const int* get_r(){ return r_; }

	const MatrixBlock* get_jaco_d() const{ return jaco_d_; }
	const MatrixBlock* get_jaco_ld() const{ return jaco_ld_; }
	const MatrixBlock* get_jaco_ud() const{ return jaco_ud_; }

	const double* get_resid() const{ return resid_; }

public:
	// Functions for test phase 3a
	void InputData();
	void OutputJacobian();
	void GenerateRestartFile();
private:
    //==========================================================
	int n_act_cell_ , nconns_;
//	int num_comps_;
	vector<string> comps_name;
	int num_matrix_blocks_;
	int hasOil, hasWater, hasGas;
	int num_phase_;
    SimCtrl::PHASETYPE phase_type_;
    //==========================================================
    //Block-wise properties
    //InitialState

	//State Variables 
    double *po_, *sw_, *sg_; //Primary Variable at n+1 time step
	double *pw_, *pg_, *so_;

    double *so_n_, *sw_n_, *sg_n_; //Main Variable at n time step
    double *pw_n_, *pg_n_, *po_n_;

    // Non-derivative properties
    double *poro_, *poro_n_;
	double *bo_, *bg_, *bw_;
    double *bo_n_,*bg_n_,*bw_n_;
    double *kro_, *krg_, *krw_;
	double *vo_, *vg_, *vw_;
    double *deno_, *deng_, *denw_;
	double deno_std_, denw_std_, deng_std_;
	double* rso_, *rso_n_;

    // Derivative properties
    double *dbo_dpo_, *dbg_dpg_, *dbw_dpw_;
    double *dkro_dsg_, *dkro_dsw_, *dkrw_dsw_, *dkrg_dsg_;
    double *dvo_dpo_, *dvg_dpg_, *dvw_dpw_;
    double *drso_dpo_;
	double *dporo_dpo_;

    // Aggregated non-derivative properties
    double *mo_,*mw_,*mg_;
    double *ao_, *aw_, *ag_;

	// Aggregated derivative properties
	double *dmo_dpo_, *dmg_dpg_, *dmw_dpw_;
	double *dmo_dsw_, *dmo_dsg_, *dmg_dsg_, *dmw_dsw_;
	double *dgammao_dpo_, *dgammag_dpg_, *dgammaw_dpw_;

    //=========================================================
    //Connection-wise properties
    const int* l_;
    const int* r_;
    const double* trans_;
    // Non-derivative property
    double *dphio_,*dphiw_,*dphig_;
    double *gammao_,*gammaw_,*gammag_;
    double *fo_, *fw_, *fg_; // Flux between connection

    //=========================================================
    // Newton-raphson variables
    // Residual vector
    double *resid_;
    //Order is Residual is Water-Oil-Gas
    // Jacobian Matrix
    // Order is variable is P,Sw,Sg
    MatrixBlock* jaco_d_; // Block-wise
    MatrixBlock* jaco_ud_; // Connection-wise
    MatrixBlock* jaco_ld_; // Connection-wise

    H5File *state_report_file_;
    DataSpace *state_report_data_space_;
    int state_report_;


};

#endif

