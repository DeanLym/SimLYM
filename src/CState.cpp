
#include"CState.h"
#include"CGrid.h"
#include"CPVT.h"
#include"CSAT.h"
//#include"CSimCtrl.h"
#include <math.h>
#include "PWLI.h"
#include "CSchedule.h"
#include <fstream>
#include <string>
#include <sstream>
#include "CStandardWell.h"
using namespace std;

//using::stringstream;


CState::CState(SimCtrl::PHASETYPE phase, CGrid *Grid, CPVT* PVT){
	n_act_cell_ = Grid->get_ncell();
	nconns_ = Grid->get_nconns();
	const ConnList *conn_temp = Grid->GetConnList();
	l_ = conn_temp->get_l();
	r_ = conn_temp->get_r();
	trans_ = conn_temp->get_trans();

	deno_std_ = PVT->GetDenO_Std();
	deng_std_ = PVT->GetDenG_Std();
	denw_std_ = PVT->GetDenW_Std();
	//	num_comps_ = num_comps;

    state_report_file_ = NULL;
    state_report_ = 1;

	phase_type_ = phase;

	switch (phase)
	{
	case SimCtrl::OWG:
		hasOil = 1;
		hasWater = 1;
		hasGas = 1;
		break;
	case SimCtrl::OW:
		hasOil = 1;
		hasWater = 1;
		hasGas = 0;
		break;
	case SimCtrl::OG:
		hasOil = 1;
		hasWater = 0;
		hasGas = 1;
		break;
	default:
		break;
	}

	num_phase_ = hasOil + hasWater + hasGas;
	resid_ = new double[num_phase_*n_act_cell_];
	memset(resid_, 0, num_phase_*n_act_cell_*sizeof(double));
	jaco_d_ = new MatrixBlock[n_act_cell_];
	jaco_ud_ = new MatrixBlock[nconns_];
	jaco_ld_ = new MatrixBlock[nconns_];

	MatrixBlock* temp;
	for(int i =0;i<n_act_cell_;i++){
		temp = new MatrixBlock(num_phase_,num_phase_);
		jaco_d_[i] = *temp;
	}
	for(int i =0;i<nconns_;i++){
		temp = new MatrixBlock(num_phase_,num_phase_);
		jaco_ud_[i] = *temp;
	}
	for(int i=0;i<nconns_;i++){
		temp = new MatrixBlock(num_phase_,num_phase_);
		jaco_ld_[i] = *temp;
	}


	poro_ = new double[n_act_cell_];
	memset(poro_, 0, n_act_cell_*sizeof(double));

	poro_n_ = new double[n_act_cell_];
	memset(poro_n_, 0, n_act_cell_*sizeof(double));

	dporo_dpo_ = new double[n_act_cell_];
	memset(dporo_dpo_, 0, n_act_cell_*sizeof(double));

	if (hasOil == 1){
		bo_ = new double[n_act_cell_];
		bo_n_ = new double[n_act_cell_];
		kro_ = new double[n_act_cell_];
		vo_ = new double[n_act_cell_];
		po_ = new double[n_act_cell_];
		po_n_ = new double[n_act_cell_];
		so_ = new double[n_act_cell_];
		so_n_ = new double[n_act_cell_];
		mo_ = new double[n_act_cell_];
		deno_ = new double[n_act_cell_];
		ao_ = new double[n_act_cell_];


		dbo_dpo_ = new double[n_act_cell_];
		dvo_dpo_ = new double[n_act_cell_];

		dmo_dpo_ = new double[n_act_cell_];
		dgammao_dpo_ = new double[n_act_cell_];


		gammao_ = new double[nconns_];
		dphio_ = new double[nconns_];
		fo_ = new double[nconns_];
	}
	if (hasGas == 1){
		bg_ = new double[n_act_cell_];
		bg_n_ = new double[n_act_cell_];
		krg_ = new double[n_act_cell_];
		vg_ = new double[n_act_cell_];
		pg_ = new double[n_act_cell_];
		pg_n_ = new double[n_act_cell_];
		sg_ = new double[n_act_cell_];
		sg_n_ = new double[n_act_cell_];
		mg_ = new double[n_act_cell_];
		deng_ = new double[n_act_cell_];
		ag_ = new double[n_act_cell_];

		dbg_dpg_ = new double[n_act_cell_];
		dvg_dpg_ = new double[n_act_cell_];
		dkrg_dsg_ = new double[n_act_cell_];


		dmg_dpg_ = new double[n_act_cell_];
		dmg_dsg_ = new double[n_act_cell_];
		dgammag_dpg_ = new double[n_act_cell_];

		gammag_ = new double[nconns_];
		dphig_ = new double[nconns_];
		fg_ = new double[nconns_];

	}
	if (hasOil == 1 && hasGas == 1){
		rso_ = new double[n_act_cell_];
		rso_n_ = new double[n_act_cell_];

		drso_dpo_ = new double[n_act_cell_];
		dkro_dsg_ = new double[n_act_cell_];

		dmo_dsg_ = new double[n_act_cell_];
	}
	if (hasOil == 1 && hasWater == 1){
		dkro_dsw_ = new double[n_act_cell_];
		dmo_dsw_ = new double[n_act_cell_];
	}
	if (hasWater == 1){
		bw_ = new double[n_act_cell_];
		bw_n_ = new double[n_act_cell_];
		krw_ = new double[n_act_cell_];
		vw_ = new double[n_act_cell_];
		pw_ = new double[n_act_cell_];
		pw_n_ = new double[n_act_cell_];
		sw_ = new double[n_act_cell_];
		sw_n_ = new double[n_act_cell_];
		mw_ = new double[n_act_cell_];
		denw_ = new double[n_act_cell_];
		aw_ = new double[n_act_cell_];

		dbw_dpw_ = new double[n_act_cell_];
		dvw_dpw_ = new double[n_act_cell_];
		dkrw_dsw_ = new double[n_act_cell_];

		dmw_dpw_ = new double[n_act_cell_];
		dmw_dsw_ = new double[n_act_cell_];
		dgammaw_dpw_ = new double[n_act_cell_];

		gammaw_ = new double[nconns_];
		dphiw_ = new double[nconns_];
		fw_ = new double[nconns_];

	}
}

CState::~CState(){
    delete state_report_file_;
    delete state_report_data_space_;

	delete[] resid_;
	delete[] jaco_d_;
	delete[] jaco_ud_;
	delete[] jaco_ld_;


	if (hasOil == 1){
		delete[] bo_;
		delete[] bo_n_;
		delete[] kro_;
		delete[] vo_;
		delete[] po_;
		delete[] po_n_;
		delete[] so_;
		delete[] so_n_;
		delete[] mo_;
		delete[] deno_;
		delete[] ao_;

		delete[] dbo_dpo_;
		delete[] dvo_dpo_;

		delete[] dmo_dpo_;
		delete[] dgammao_dpo_;

		delete[] gammao_;
		delete[] dphio_;
		delete[] fo_;

	}
	if (hasGas == 1){

		delete[] bg_;
		delete[] bg_n_;
		delete[] krg_;
		delete[] vg_;
		delete[] pg_;
		delete[] pg_n_;
		delete[] sg_;
		delete[] sg_n_;
		delete[] mg_;
		delete[] deng_;
		delete[] ag_;

		delete[] dbg_dpg_;
		delete[] dvg_dpg_;
		delete[] dkrg_dsg_;

		delete[] dmg_dpg_;
		delete[] dmg_dsg_;
		delete[] dgammag_dpg_;

		delete[] gammag_;
		delete[] dphig_;
		delete[] fg_;
	}
	if (hasOil == 1 && hasGas == 1){

		delete[] rso_;
		delete[] rso_n_;

		delete[] drso_dpo_;
		delete[] dkro_dsg_;
		delete[] dmo_dsg_;

	}
	if (hasOil == 1 && hasWater == 1){
		delete[] dkro_dsw_;
		delete[] dmo_dsw_;
	}
	if (hasWater == 1){
		delete[] bw_;
		delete[] krw_;
		delete[] vw_;
		delete[] pw_;
		delete[] sw_;
		delete[] mw_;
		delete[] denw_;
		delete[] aw_;

		delete[] dbw_dpw_;
		delete[] dvw_dpw_;
		delete[] dkrw_dsw_;

		delete[] dmw_dpw_;
		delete[] dmw_dsw_;
		delete[] dgammaw_dpw_;

		delete[] gammaw_;
		delete[] dphiw_;
		delete[] fw_;
	}
	delete[] poro_;
	delete[] poro_n_;
	delete[] dporo_dpo_;

}

bool CState::SetInitPres(double po){
	for (int m = 0; m < n_act_cell_; m++){
		po_[m] = pw_[m] = po;
	}

	return true;
}

bool CState::SetInitPres(double *po){
	for (int m = 0; m < n_act_cell_; m++){
		po_[m] = pw_[m] = po[m];
	}
	return true;
}

bool CState::SetInitPres(char *po_file){
	ifstream in;
	in.open(po_file);
	if (in.is_open()){
		for (int i = 0; i < n_act_cell_; i++){
			in >> po_[i];
			pw_[i] = po_[i];
		}
	}
	else{
		cout << "Open Initial Pressure File Failed!" << endl;
		exit(1);
		return false;
	}
	in.close();
	return true;
}

bool CState::SetInitPres(double datumDepth, double pDatum, CPVT *PVT, CGrid *Grid){

	//Construct a pressure table for interpolation;
	double dD;
	dD = 6; //6 ft between each pressure table point;
	int nNodeDown, nNodeUp;

	nNodeDown = ceil((Grid->GetdMax() - datumDepth) / dD) + 1;
	nNodeUp = ceil((datumDepth - Grid->GetdMin()) / dD) + 1;
	//Downward

	double *dDown, *dUp;
	double *pDown, *pUp;
	double pOld, pNew, denMean;
	double err = 1e-6;
	// i = 0;
	if (nNodeDown>1){
		dDown = new double[nNodeDown];
		for (int i = 0; i < nNodeDown; i++){
			dDown[i] = datumDepth + dD*i;
		}
		pOld = pDatum + PVT->GetDenO(pDatum)*dD*beta;
		denMean = (PVT->GetDenO(pDatum) + PVT->GetDenO(pOld)) / 2;
		pNew = pDatum + denMean*dD*beta;

		while (fabs(pOld - pNew) > err){
			pOld = pNew;
			denMean = PVT->GetDenO((pDatum + pNew) / 2);
			pNew = pDatum + denMean*dD*beta;
		}
		pDown = new double[nNodeDown];
		pDown[0] = pDatum;
		pDown[1] = pNew;
		for (int i = 2; i < nNodeDown; i++){
			pOld = pDown[i - 1] + PVT->GetDenO(pDown[i - 1])*dD*beta;
			denMean = (PVT->GetDenO(pDown[i - 1]) + PVT->GetDenO(pOld)) / 2;
			pNew = pDown[i - 1] + denMean*dD*beta;
			while (fabs(pOld - pNew) > err){
				pOld = pNew;
				denMean = (PVT->GetDenO(pDown[i - 1]) + PVT->GetDenO(pNew)) / 2;
				pNew = pDown[i - 1] + denMean*dD*beta;
			}
			pDown[i] = pNew;
		}
	}
	else{
		dDown = new double[2];
		dDown[0] = datumDepth;
		dDown[1] = datumDepth;
		pDown = new double[2];
		pDown[0] = pDatum;
		pDown[1] = pDatum;
		nNodeDown = 2;
	}

	//Upward
	if (nNodeUp > 1){
		dUp = new double[nNodeUp];
		for (int i = 0; i < nNodeUp; i++){
			dUp[i] = datumDepth - dD*i;
		}
		pOld = pDatum - PVT->GetDenO(pDatum)*dD*beta;
		denMean = (PVT->GetDenO(pDatum) + PVT->GetDenO(pOld)) / 2;
		pNew = pDatum - denMean*dD*beta;
		while (fabs(pOld- pNew) > err){
			pOld = pNew;
			denMean = (PVT->GetDenO(pDatum) + PVT->GetDenO(pNew)) / 2;
			pNew = pDatum - denMean*dD*beta;
		}
		pUp = new double[nNodeUp];
		pUp[0] = pDatum;
		pUp[1] = pNew;
		for (int i = 1; i < nNodeUp; i++){
			pOld = pUp[i - 1] - PVT->GetDenO(pUp[i - 1])*dD*beta;
			denMean = (PVT->GetDenO(pUp[i - 1]) + PVT->GetDenO(pOld)) / 2;
			pNew = pUp[i - 1] - denMean*dD*beta;
			while (fabs(pOld - pNew) > err){
				pOld = pNew;
				denMean = (PVT->GetDenO(pUp[i - 1]) + PVT->GetDenO(pNew)) / 2;
				pNew = pUp[i - 1] + denMean*dD*beta;
			}
			pUp[i] = pNew;
		}
	}
	else{
		dUp = new double[2];
		dUp[0] = datumDepth;
		dUp[1] = datumDepth;
		pUp = new double[2];
		pUp[0] = pDatum;
		pUp[1] = pDatum;
		nNodeUp = 2;
	}


	int m;
	double depM;
	for (int m = 0; m < n_act_cell_; m++){
		depM = Grid->GetDepth(m);
		if (depM>datumDepth){
			po_[m] = PiecewiseLinearInter(nNodeDown, dDown, pDown, depM);
		}
		else{
			if (depM < datumDepth){
				po_[m] = PiecewiseLinearInter(nNodeUp, dUp, pUp, depM);
			}
			else{
				po_[m] = pDatum;
			}
		}
		pw_[m] = po_[m];
		//	pg_[m] = pg_n_[m] = po_[m];
	}
	delete[] pDown;
	delete[] pUp;
	delete[] dDown;
	delete[] dUp;
	//while 
	return 1;

}

bool CState::SetInitSat(double Sw, double Sg){
	// Only for simple cases with homogenous initial saturation field.
	double So;
	So = 1 - Sw - Sg;
	for (int i = 0; i < n_act_cell_; i++){
		sw_[i] = Sw;
		//		sg_[i] = sg_n_[i] = Sg;
		so_[i] = So;
	}
	return 1;
}

bool CState::SetInitSat(double *Sw, double *Sg){
	// Only for simple cases with homogenous initial saturation field.
	double So;

	for (int i = 0; i < n_act_cell_; i++){
		So = 1 - Sw[i] - Sg[i];
		sw_[i] = Sw[i];
		so_[i] = So;
	}
	return 1;
}

bool CState::SetInitSat(char *Sw_file, char *Sg_file){
	ifstream in;
	in.open(Sw_file);
	if (in.is_open()){
		for (int i = 0; i < n_act_cell_; i++){
			in >> sw_[i];
		}
	}
	else{
		cout << "Open Initial Sw File Failed!" << endl;
		exit(1);
		return false;
	}
	in.close();

	in.open(Sg_file);
	if (in.is_open()){
		for (int i = 0; i < n_act_cell_; i++){
			in >> sg_[i];
			so_[i] = 1 - sw_[i] - sg_[i];
		}
	}
	else{
		cout << "Open Initial Sg File Failed!" << endl;
		exit(1);
		return false;
	}
	in.close();
	return true;
}


bool CState::SetInitSw(char *Sw_file){
	ifstream in;
	in.open(Sw_file);
	if (in.is_open()){
		for (int i = 0; i < n_act_cell_; i++){
			in >> sw_[i];
			so_[i] = 1 - sw_[i];
		}
	}
	else{
		cout << "Open Initial Sw File Failed!" << endl;
		exit(1);
		return false;
	}
	in.close();
	return true;
}

bool CState::SetInitSg(char *Sg_file){
	ifstream in;
	in.open(Sg_file);
	if (in.is_open()){
		for (int i = 0; i < n_act_cell_; i++){
			in >> sg_[i];
			so_[i] = 1 - sg_[i];
		}
	}
	else{
		cout << "Open Initial Sw File Failed!" << endl;
		exit(1);
		return false;
	}
	in.close();
	return true;
}


bool CState::CalDPhiOil(CGrid* Grid){
	// Calculate Potential Difference;
	const double* delta_d_ = Grid->GetDeltaD();
	for (int i = 0; i < nconns_; i++)
		dphio_[i] = po_[l_[i]] - po_[r_[i]] - gammao_[i] * delta_d_[i];
	return true;
}

bool CState::CalDPhiGas(CGrid* Grid){
	// Calculate Potential Difference;
	const double* delta_d_ = Grid->GetDeltaD();
	for (int i = 0; i < nconns_; i++)
		dphig_[i] = pg_[l_[i]] - pg_[r_[i]] - gammag_[i] * delta_d_[i];
	return true;
}

bool CState::CalDPhiWater(CGrid* Grid){
	// Calculate Potential Difference;
	const double* delta_d_ = Grid->GetDeltaD();
	for (int i = 0; i < nconns_; i++)
		dphiw_[i] = pw_[l_[i]] - pw_[r_[i]] - gammaw_[i] * delta_d_[i];
	return true;
}

bool CState::CalGammaOil(){
	for (int i = 0; i < nconns_; i++)
		gammao_[i] = beta*(g_/gc)*(deno_[l_[i]]+deno_[r_[i]])/2;
	return true;
}

bool CState::CalGammaGas(){
	for (int i = 0; i < nconns_; i++)
		gammag_[i] = beta*(g_ / gc)*(deng_[l_[i]] + deng_[r_[i]]) / 2;
	return true;
}

bool CState::CalGammaWater(){
	for (int i = 0; i < nconns_; i++)
		gammaw_[i] = beta*(g_ / gc)*(denw_[l_[i]] + denw_[r_[i]]) / 2;
	return true;
}

bool CState::CalMobiOil(){
	int temp;
	for(int i=0;i<n_act_cell_;i++){
		mo_[i] = kro_[i]/(vo_[i]*bo_[i]);
	}
	return true;
}

bool CState::CalMobiGas(){
	for(int i=0;i<n_act_cell_;i++){
		mg_[i] = krg_[i]/(vg_[i]*bg_[i]);
	}
	return true;
}

bool CState::CalMobiWater(){
	for(int i=0;i<n_act_cell_;i++){
		mw_[i] = krw_[i]/(vw_[i]*bw_[i]);
	}
	return true;
}

bool CState::CalAccuOil(CGrid *Grid, double dt){
	const double* v;
	v = Grid->GetVolume();
	Grid->GetPoro(po_, poro_);
	Grid->GetPoro(po_n_, poro_n_);
	for (int i = 0; i < n_act_cell_; i++){
		ao_[i] = (1 / M)*v[i]*((poro_[i]* so_[i] / bo_[i]) - (poro_n_[i] * so_n_[i] / bo_n_[i])) / dt;
	}
	return true;
}

bool CState::CalAccuWater(CGrid *Grid, double dt){
	const double* v;
	v = Grid->GetVolume();
	Grid->GetPoro(po_, poro_);
	Grid->GetPoro(po_n_, poro_n_);
	for (int i = 0; i < n_act_cell_; i++){
		aw_[i] = (1 / M)*v[i]*((poro_[i] * sw_[i] / bw_[i]) - (poro_n_[i] * sw_n_[i] / bw_n_[i])) / dt;
	}
	return true;
}

bool CState::CalAccuGas(CGrid* Grid, double dt){
	const double* v;
	v = Grid->GetVolume();
	Grid->GetPoro(po_, poro_);
	Grid->GetPoro(po_n_, poro_n_);
	for (int i = 0; i < n_act_cell_; i++){
		ag_[i] = v[i]/ (M*dt)
            		*(sg_[i] / bg_[i]*poro_[i]  - sg_n_[i] / bg_n_[i] * poro_n_[i]
            		                                     + so_[i] * rso_[i] / bo_[i] * poro_[i]  - so_n_[i] * rso_n_[i] / bo_n_[i] * poro_n_[i]);
	}
	return true;
}

void CState::UpdateNewState(const double* dx){
	if (hasOil == 1 && hasGas == 1 && hasWater == 0){
		for (int n = 0; n < n_act_cell_; n++)
		{
			po_[n] += dx[2 * n];
			pg_[n] = po_[n];

		}
		for (int n = 0; n < n_act_cell_; n++)
		{
			sg_[n] += dx[2 * n + 1];
			if (sg_[n] < 0){
				sg_[n] = 0;
			}
			else if (sg_[n]>1){
				sg_[n] = 1;
			}
			so_[n] = 1 - sg_[n];
		}
	}
	else if (hasOil == 1 && hasGas == 0 && hasWater == 1){
		for (int n = 0; n < n_act_cell_; n++)
		{
			po_[n] += dx[2 * n];
			pw_[n] = po_[n];
		}
		for (int n = 0; n < n_act_cell_; n++)
		{
			sw_[n] += dx[2 * n + 1];
			if (sw_[n] < 0){
				sw_[n] = 0;
			}
			else if (sw_[n]>1){
				sw_[n] = 1;
			}
			so_[n] = 1 - sw_[n];
		}
	}
}

bool CState::UpdateOldState(CPVT *PVT){
	if (hasOil == 1 && hasGas == 1 && hasWater == 0){
		for (int n = 0; n < n_act_cell_; n++)
		{
			po_n_[n] = po_[n];
			pg_n_[n] = po_n_[n];
		}
		for (int n = 0; n < n_act_cell_; n++)
		{
			sg_n_[n] = sg_[n];
			so_n_[n] = 1 - sg_n_[n];
		}
		PVT->GetBo(n_act_cell_,po_n_,bo_n_);
		PVT->GetBg(n_act_cell_,pg_n_,bg_n_);
		PVT->GetRso(n_act_cell_,po_n_,rso_n_);
	}
	else if (hasOil == 1 && hasGas == 0 && hasWater == 1){
		for (int n = 0; n < n_act_cell_; n++)
		{
			po_n_[n] = po_[n];
			pw_n_[n] = po_n_[n];
		}
		for (int n = 0; n < n_act_cell_; n++)
		{
			sw_n_[n] = sw_[n];
			so_n_[n] = 1 - sw_n_[n];
		}
		PVT->GetBo(n_act_cell_, po_n_, bo_n_);
		PVT->GetBw(n_act_cell_, pw_n_, bw_n_);
	}
	return true;
}

void CState::ChangeBackState(){//When convergence fails
	if (hasOil == 1 && hasGas == 1 && hasWater == 0){
		for (int n = 0; n < n_act_cell_; n++)
		{
			po_[n] = po_n_[n];
			pg_[n] = po_[n];
		}
		for (int n = 0; n < n_act_cell_; n++)
		{
			sg_[n] = sg_n_[n];
			so_[n] = 1 - sg_[n];
		}
	}if (hasOil == 1 && hasGas == 0 && hasWater == 1){
		for (int n = 0; n < n_act_cell_; n++)
		{
			po_[n] = po_n_[n];
			pw_[n] = po_[n];
		}
		for (int n = 0; n < n_act_cell_; n++)
		{
			sw_[n] = sw_n_[n];
			so_[n] = 1 - sw_[n];
		}

	}
}

//void CState::State_Report(CSchedule *SCH, int n){
//	ostringstream abc;
//	abc << n;
//	string file_name;
//	file_name = abc.str();
//	ofstream out(file_name.c_str());
//	out << "Time:" << endl;
//	out << SCH->GetTCurrent() << endl;
//	out << "Po:" << endl;
//	for (int n = 0; n < n_act_cell_; n++){
//		out << n << "  " << po_[n] << endl;
//	}
//	out << "Sw:" << endl;
//	for (int n = 0; n < n_act_cell_; n++){
//		out << n << "  " << sw_[n] << endl;
//	}
//	out << "Kro:" << endl;
//	for (int n = 0; n < n_act_cell_; n++){
//		out << n << "  " << kro_[n] << endl;
//	}
//	out << "Krw:" << endl;
//	for (int n = 0; n < n_act_cell_; n++){
//		out << n << "  " << krw_[n] << endl;
//	}
//	out << "Bo:" << endl;
//	for (int n = 0; n < n_act_cell_; n++){
//		out << n << "  " << bo_[n] << endl;
//	}
//	out << "Bw:" << endl;
//	for (int n = 0; n < n_act_cell_; n++){
//		out << n << "  " << bw_[n] << endl;
//	}
//	out << "Vo:" << endl;
//	for (int n = 0; n < n_act_cell_; n++){
//		out << n << "  " << vo_[n] << endl;
//	}
//	out << "Vw:" << endl;
//	for (int n = 0; n < n_act_cell_; n++){
//		out << n << "  " << vw_[n] << endl;
//	}
//	out << "mw:" << endl;
//	for (int n = 0; n < n_act_cell_; n++){
//		out << n << "  " << mw_[n] << endl;
//	}
//	out << "mo:" << endl;
//	for (int n = 0; n < n_act_cell_; n++){
//		out << n << "  " << mo_[n] << endl;
//	}
//
//
//	/*
//	out << "Rso:" << endl;
//	for (int n = 0; n < n_act_cell_; n++){
//
//		out << n << "  " << rso_[n] / 178.1 << endl;
//	}*/
//}

void CState::Init_State_Report(int state_report){
    state_report_ = state_report;
    if (state_report == 3 || state_report_ == 2){
        state_report_file_ = new H5File("state.h5", H5F_ACC_TRUNC);
        const int RANK = 1;
        hsize_t dims[1];
        dims[0] = n_act_cell_;
        state_report_data_space_ = new DataSpace(RANK, dims);
    }
}

void CState::State_Report(CSchedule *SCH, int n){
    ostringstream abc;
    abc << "p";
    abc << n;
    string file_name;
    file_name = abc.str();
    ostringstream ab;
    ab << "sw";
    ab << n;
    string file_name2;
    file_name2 = ab.str();
    // cout << state_report_ << endl;
    if(state_report_ == 3 || state_report_ == 2){
        // cout << "Writing HDF5 file" << endl;

        const H5std_string ds_name_p(file_name.c_str());
        DataSet dataset_p = state_report_file_->createDataSet(ds_name_p, PredType::NATIVE_DOUBLE, *state_report_data_space_);
        dataset_p.write(po_, PredType::NATIVE_DOUBLE);

        const H5std_string ds_name_sw(file_name2.c_str());
        DataSet dataset_sw = state_report_file_->createDataSet(ds_name_sw, PredType::NATIVE_DOUBLE, *state_report_data_space_);
        dataset_sw.write(sw_, PredType::NATIVE_DOUBLE);

    }else{
        
        ofstream out(file_name.c_str());
        for (int n = 0; n < n_act_cell_; n++){
            out << po_[n] << endl;
        }
        out.close();


        ofstream out2(file_name2.c_str());
        for (int n = 0; n < n_act_cell_; n++){
            out2 << sw_[n] << endl;
        }
        out2.close();        
    }


}

void CState::State_Debug(){
	ofstream out("state_debug.txt");
	out << "Po:" << endl;
	for (int n = 0; n < n_act_cell_; n++){
		out << n << "  " << po_[n] << endl;
	}
	out << "So:" << endl;
	for (int n = 0; n < n_act_cell_; n++){
		out << n << "  " << so_[n] << endl;
	}
	/*
	out << "Rso:" << endl;
	for (int n = 0; n < n_act_cell_; n++){

		out << n << "  " << rso_[n] / 178.1 << endl;
	}
	 */
}


void CState::GenerateRestartFile(){
	if(hasOil){
		ofstream out1("po.txt");
		for (int n = 0; n < n_act_cell_; n++){
			out1 << po_[n] << endl;
		}
		out1.close();
	}
	if(hasWater){
		ofstream out2("sw.txt");
		for (int n = 0; n < n_act_cell_; n++){
			out2 << sw_[n] << endl;
		}
		out2.close();
	}
	if(hasGas){
		ofstream out3("sg.txt");
		for (int n = 0; n < n_act_cell_; n++){
			out3 << sg_[n] << endl;
		}
		out3.close();
	}

}

void CState::Var_Debug(double *var, int n, char* fn){
	ofstream out(fn);
	for (int i = 0; i < n; i++){
		out << i << "  " << var[i] << endl;
	}
}

void CState::Residual_Debug(char* fn){
	ofstream out(fn);
	out << "Residual:" << endl;
	for (int n = 0; n < num_phase_ * n_act_cell_; n++){
		out << n << "  " << resid_[n] << endl;
	}
}

bool CState::CalFo(){
	double mobio = 0.0;

	for(int i = 0;i<nconns_;i++){
		if(dphio_[i] > 0){
			mobio = mo_[l_[i]];
		}else{
			mobio = mo_[r_[i]];
		}
		fo_[i] = trans_[i]*mobio*dphio_[i];
	}
	return true;
}

bool CState::AssembleResidualOG(vector<CStandardWell*> &std_well){
	//Assemble Residual Vector
	int index_r = 0;
	int index_l = 0;
	int n = 0;
	for (int i = 0; i < n_act_cell_; i++){
		resid_[2 * i] = -ao_[i];
		resid_[2 * i + 1] = -ag_[i];
	}
	for (int i = 0; i<nconns_; i++){
		// residual is Ro, Rg
		index_r = 2 * r_[i];
		index_l = 2 * l_[i];
		resid_[index_r] += fo_[i];
		resid_[index_l] -= fo_[i];
		resid_[index_r + 1] += fg_[i];
		resid_[index_l + 1] -= fg_[i];
	}
	for (int i = 0; i < std_well.size(); i++){
		n = std_well[i]->get_n();
		resid_[2 * n] -= std_well[i]->get_qo();
		resid_[2 * n + 1] -= std_well[i]->get_qg();
	}
	return true;
}

bool CState::AssembleResidualOW(vector<CStandardWell*> &std_well){
	int index_r = 0;
	int index_l = 0;
	int n = 0;
	for (int i = 0; i < n_act_cell_; i++){
		resid_[2 * i] = -aw_[i]; // Add up accumulation term of water to Rw(i)
		resid_[2 * i + 1] = -ao_[i]; //Add up accumulation term of oil to Ro(i)
	}
	//	Residual_Debug("residual_accumulation.txt");
	for(int i = 0;i<nconns_;i++){
		// residual is Rw, Ro
		index_r = 2*r_[i]; // Grid block index of the right side of the connection
		index_l = 2*l_[i]; // Grid block index of the left side of the connection
		resid_[index_r] += fw_[i]; // Add up flux term of water to Rw(r)
		resid_[index_l] -= fw_[i]; // Add up flux term of water to Rw(l)
		resid_[index_r+1] += fo_[i]; // Add up flux  term of oil to Ro(r)
		resid_[index_l+1] -= fo_[i]; // Add up flux term of oil to Ro(l)
	}
	//	Residual_Debug("residual_accum_flux.txt");
	for (int i = 0; i < std_well.size(); i++){
		n = std_well[i]->get_n(); //cout << "n: "<< n ;
		resid_[2 * n] -= std_well[i]->get_qw(); //cout<<"qw:"<<std_well[i]->get_qw();
		resid_[2 * n + 1] -= std_well[i]->get_qo(); //cout<<"qo:"<<std_well[i]->get_qo();
	}
	//	Residual_Debug("residual_all.txt");
	//	system("pause");
	return true;
}

bool CState::AssembleResidualOWG(vector<CStandardWell*> &std_well){
	int index_r = 0;
	int index_l = 0;
	int n = 0;
	for (int i = 0; i < n_act_cell_; i++){
		resid_[3 * i] = -aw_[i];
		resid_[3 * i + 1] = -ao_[i];
		resid_[3 * i + 2] = -ag_[i];
	}
	for(int i = 0;i<nconns_;i++){
		// residual is Rw, Ro, Rg
		index_r = 3*r_[i];
		index_l = 3*l_[i];
		resid_[index_r] += fw_[i];
		resid_[index_l] -= fw_[i];
		resid_[index_r+1] += fo_[i];
		resid_[index_l+1] -= fo_[i];
		resid_[index_r+2] += fg_[i];
		resid_[index_l+2] -= fg_[i];
	}
	for (int i = 0; i < std_well.size(); i++){
		n = std_well[i]->get_n();
		resid_[3 * n] -= std_well[i]->get_qw();
		resid_[3 * n + 1] -= std_well[i]->get_qo();
		resid_[3 * n + 2] -= std_well[i]->get_qg();
	}
	return true;
}

bool CState::CalFw(){
	double mobiw = 0.0;
	for(int i = 0;i<nconns_;i++){
		if(dphiw_[i] > 0){
			mobiw = mw_[l_[i]];
		}else{
			mobiw = mw_[r_[i]];
		}
		//        if(dphiw_ >)
			fw_[i] = trans_[i]*mobiw*dphiw_[i];
	}
	return true;
}

bool CState::CalFg(){
	double mobig = 0.0;
	double rso = 0.0;
	for(int i = 0;i<nconns_;i++){
		if(dphig_[i] > 0){
			mobig = mg_[l_[i]];
		}else{
			mobig = mg_[r_[i]];
		}
		if(dphio_[i] > 0){
			rso = rso_[l_[i]];
		}else{
			rso = rso_[r_[i]];
		}
		fg_[i] = trans_[i]*(mobig*dphig_[i]+rso*fo_[i]);
	}
	return true;
}

bool CState::AssembleJacobianOG(CGrid* Grid, double dt, vector<CStandardWell*> &std_well){
	MatrixBlock* temp;
	for (int i = 0; i<n_act_cell_; i++){

	}
	return true;
}

bool CState::CalDMwDPw(){
	for (int i = 0; i < n_act_cell_; i++){
		dmw_dpw_[i] = -krw_[i] * (dbw_dpw_[i] * vw_[i] + bw_[i] * dvw_dpw_[i]) / (vw_[i] * vw_[i] * bw_[i] * bw_[i]); 
	}
	return true;
}

bool CState::CalDMoDPo(){
	for (int i = 0; i < n_act_cell_; i++){
		dmo_dpo_[i] = -kro_[i] * (dbo_dpo_[i] * vo_[i] + bo_[i] * dvo_dpo_[i]) / (vo_[i] * vo_[i] * bo_[i] * bo_[i]);
	}
	return true;
}

bool CState::CalDMwDSw(){
	for (int i = 0; i < n_act_cell_; i++){
		dmw_dsw_[i] = dkrw_dsw_[i] / (vw_[i] * bw_[i]);
	}
	return true;
}

bool CState::CalDMoDSw(){
	for (int i = 0; i < n_act_cell_; i++){
		dmo_dsw_[i] = dkro_dsw_[i] / (vo_[i] * bo_[i]);
	}
	return true;
}

bool CState::CalDGammaODPo(){
	for (int i = 0; i < n_act_cell_; i++){
		dgammao_dpo_[i] = -0.5*beta*(g_ / gc)*deno_std_*dbo_dpo_[i] / (bo_[i] * bo_[i]);
	}
	return true;
}

bool CState::CalDGammaWDPw(){
	for (int i = 0; i < n_act_cell_; i++){
		dgammaw_dpw_[i] = -0.5*beta*(g_ / gc)*denw_std_*dbw_dpw_[i] / (bw_[i] * bw_[i]);
	}
	return true;
}

bool CState::CalDPoroDPo(CGrid *Grid){
	return Grid->GetDPoroDPo(po_, dporo_dpo_);
}

bool CState::AssembleJacobianOW(CGrid* Grid, double dt, vector<CStandardWell*> &std_well){
	MatrixBlock* temp;
	const double* v; //Pointer of volume array
	v = Grid->GetVolume();
	const double* delta_d = Grid->GetDeltaD(); //Pointer to delta_d array


	//Flux Term
	double m_lr;
	double dmw_dpl, dmw_dpr;
	double dmw_dswl, dmw_dswr;
	double dmo_dpl, dmo_dpr;
	double dmo_dswl, dmo_dswr;
	int l, r;
	//Accumulation Term
	for (int i = 0; i<n_act_cell_; i++){
		jaco_d_[i](0, 0) = -v[i] * sw_[i] * (dporo_dpo_[i] * bw_[i] - dbw_dpw_[i] * poro_[i]) / (M*dt*bw_[i] * bw_[i]); //dAw(i)/dPo(i)
		jaco_d_[i](0, 1) = -v[i] * poro_[i] / (M*dt*bw_[i]); //dAw(i)/dSw(i)

		jaco_d_[i](1, 0) = -v[i] * so_[i] * (dporo_dpo_[i] * bo_[i] - dbo_dpo_[i] * poro_[i]) / (M*dt*bo_[i] * bo_[i]);//dAo(i)/dPo(i)
		jaco_d_[i](1, 1) = v[i] * poro_[i] / (M*dt*bo_[i]);// dAo(i)/dSw(i)
	}

	for (int i = 0; i < nconns_; i++){
		l = l_[i];
		r = r_[i];

		// Water Phase
		// Upwinding treatment of mobility term
		if (dphiw_[i]>0){  // phiw(l) > phiw(r)
			m_lr = mw_[l]; //Mobility of water at connection lr, upstream weighting
			dmw_dpr = 0; //dmw(lr)/dp(r)
			dmw_dpl = dmw_dpw_[l]; //dmw(lr)/dp(l)
			dmw_dswl = dmw_dsw_[l]; //dmw(lr)/dsw(l)
			dmw_dswr = 0; //dmw(lr)/dsw(r)
		}
		else{
			m_lr = mw_[r]; //Mobility of water at connection lr, upstream weighting
			dmw_dpr = dmw_dpw_[r]; //dmw(lr)/dp(r)
			dmw_dpl = 0; //dmw(lr)/dp(l)
			dmw_dswl = 0; //dmw(lr)/dsw(l)
			dmw_dswr = dmw_dsw_[r]; //dmw(lr)/dsw(r)
		}
		jaco_d_[l](0, 0) += trans_[i] * (-dmw_dpl*dphiw_[i] + m_lr*(-1 + dgammaw_dpw_[l]*delta_d[i])); //dFw(l)/dP(l)
		jaco_d_[l](0, 1) += -trans_[i] * dmw_dswl*dphiw_[i]; //dFw(l)/dSw(l)

		jaco_d_[r](0, 0) += trans_[i] * (dmw_dpr*dphiw_[i] + m_lr*(-1 - dgammaw_dpw_[r]*delta_d[i]));  //dFw(r)/dP(r)
		jaco_d_[r](0, 1) += trans_[i] * dmw_dswr*dphiw_[i]; //dFw(r)/dSw(r)

		// l is always smaller than r by construction
		// Upper diagonal has col>row, corresponding to J(l,r)
		jaco_ud_[i](0, 0) = trans_[i] * (-dmw_dpr*dphiw_[i] + m_lr*(1 + dgammaw_dpw_[r] * delta_d[i]));  //dFw(l)/dP(r);
		jaco_ud_[i](0, 1) = -trans_[i] * dmw_dswr*dphiw_[i];//dFw(l)/dSw(r)
		// Lower diagonal has col<row, corresponding to J(r,1)
		jaco_ld_[i](0, 0) = trans_[i] * (dmw_dpl*dphiw_[i] + m_lr*(1 - dgammaw_dpw_[l] * delta_d[i]));  //dFw(r)/dP(l);
		jaco_ld_[i](0, 1) = trans_[i] * dmw_dswl*dphiw_[i];//dFw(r)/dSw(l)

		// Oil Phase
		// Upwinding treatment of mobility term
		if (dphio_[i]>0){  // phio(l) > phio(r)
			m_lr = mo_[l]; // Mobility of oil at connection between l and r
			dmo_dpr = 0; //dmo(lr)/dp(r)
			dmo_dpl = dmo_dpo_[l]; //dmo(lr)/dp(l)
			dmo_dswl = dmo_dsw_[l]; //dmo(lr)/dsw(l)
			dmo_dswr = 0; //dmo(lr)/dsw(r)
		}
		else{
			m_lr = mo_[r]; // Mobility of oil at connection between l and r
			dmo_dpr = dmo_dpo_[r];  //dmo(lr)/dp(r)
			dmo_dpl = 0; //dmo(lr)/dp(l)
			dmo_dswl = 0; //dmo(lr)/dsw(l)
			dmo_dswr = dmo_dsw_[r]; //dmo(lr)/dsw(r)
		} 
		jaco_d_[l](1, 0) += trans_[i] * (-dmo_dpl*dphio_[i] + m_lr*(-1 + dgammao_dpo_[l]*delta_d[i])); //dFo(l)/dP(l)
		jaco_d_[l](1, 1) += -trans_[i] * dmo_dswl*dphio_[i]; //dFo(l)/dSw(l)

		jaco_d_[r](1, 0) += trans_[i] * (dmo_dpr*dphio_[i] + m_lr*(-1 - dgammao_dpo_[r]*delta_d[i]));  //dFo(r)/dP(r)
		jaco_d_[r](1, 1) += trans_[i] * dmo_dswr*dphio_[i]; //dFo(r)/dSw(r)

		// l is always smaller than r by construction
		// Upper diagonal has col>row, corresponding to J(l,r)
		jaco_ud_[i](1, 0) = trans_[i] * (-dmo_dpr*dphio_[i] + m_lr*(1 + dgammao_dpo_[r] * delta_d[i]));  //dFo(l)/dP(r);
		jaco_ud_[i](1, 1) = -trans_[i] * dmo_dswr*dphio_[i];//dFo(l)/dSw(r)
		// Lower diagonal has col<row, corresponding to J(r,1)
		jaco_ld_[i](1, 0) = trans_[i] * (dmo_dpl*dphio_[i] + m_lr*(1 - dgammao_dpo_[l] * delta_d[i])); //dFo(r)/dP(l);
		jaco_ld_[i](1, 1) = trans_[i] * dmo_dswl*dphio_[i];//dFo(r)/dSw(l)
	}
	int n = 0;
	for (int i = 0; i < std_well.size(); i++){
		n = std_well[i]->get_n();
		jaco_d_[n](0, 0) -= std_well[i]->get_dqw_dp();
		jaco_d_[n](0, 1) -= std_well[i]->get_dqw_dsw();
		jaco_d_[n](1, 0) -= std_well[i]->get_dqo_dp();
		jaco_d_[n](1, 1) -= std_well[i]->get_dqo_dsw();
	}
	return true;
}

bool CState::AssembleJacobianOWG(CGrid* Grid, double dt, vector<CStandardWell*> &std_well){
	return true;
}

void CState::CalParametersJacobian(CGrid *Grid, vector<CStandardWell*> &std_well, CPVT *PVT, CSAT *SAT, CSchedule *SCH){
	if (phase_type_ == SimCtrl::OG){

	}
	else if (phase_type_ == SimCtrl::OW){
		//=========================================================
		// Calculate the derivative properties
		PVT->GetBwPD(n_act_cell_, pw_, dbw_dpw_);  //Calculate dBw/dPw at time step n+1
		PVT->GetBoPD(n_act_cell_, po_, dbo_dpo_);  //Calculate dBo/dPo at time step n+1
		PVT->GetVisOPD(n_act_cell_, po_, dvo_dpo_); //Calculate dVo/dPo at time step n+1
		PVT->GetVisWPD(n_act_cell_, pw_, dvw_dpw_); //Calculate dVw/dPw at time step n+1
		SAT->GetKroPDSw(n_act_cell_, sw_, dkro_dsw_); //Calculate dKro/dSw at time step n+1
		SAT->GetKrwPDSw(n_act_cell_, sw_, dkrw_dsw_); //Calculate dKrw/dSw at time step n+1

		CalDMoDPo(); //Calculate dMo/dPo at time step n+1
		CalDMoDSw(); //Calculate dMo/dSw at time step n+1
		CalDMwDPw(); //Calculate dMw/dPw at time step n+1
		CalDMwDSw(); //Calculate dMw/dSw at time step n+1

		CalDGammaODPo(); //Calculate dgammao/dPo at time step n+1
		CalDGammaWDPw(); //Calculate dgammaw/dPw at time step n+1

		CalDPoroDPo(Grid);
	}
	else if (phase_type_ == SimCtrl::OWG){

	}
}

void CState::CalParametersResidual(CGrid *Grid, vector<CStandardWell*> &std_well, CPVT *PVT, CSAT *SAT, CSchedule *SCH){
	double dt = SCH->GetDt();
	if (phase_type_ == SimCtrl::OG){
		//=====================================================
		// Calculate Block-wise properties
		// Calculate the non-derivative properties
		PVT->GetBo(n_act_cell_,po_,bo_);
		PVT->GetBg(n_act_cell_,pg_,bg_);
		PVT->GetVisO(n_act_cell_,po_,vo_);
		PVT->GetVisG(n_act_cell_,pg_,vg_);
		PVT->GetRso(n_act_cell_,po_,rso_);
		PVT->GetDenO(n_act_cell_,bo_,deno_);
		PVT->GetDenG(n_act_cell_,bg_,deng_);
		SAT->GetKro(n_act_cell_,so_,kro_);
		SAT->GetKrg(n_act_cell_,sg_,krg_);

		CalGammaOil();
		CalGammaGas();

		CalMobiOil();
		CalMobiGas();

		CalAccuOil(Grid, dt);
		CalAccuGas(Grid, dt);


		//=====================================================
				// Calculate Connection-wise propeties

		// Non derivative properties
		CalGammaOil();
		CalGammaGas();
		CalDPhiOil(Grid);
		CalDPhiGas(Grid);
		CalFo();
		CalFg();

	}
	else if (phase_type_ == SimCtrl::OW){
		//=====================================================
		//  Calculate Residual Vector
		// Calculate Block-wise properties
		PVT->GetBo(n_act_cell_, po_, bo_);	//Calculate Bo at time step n+1
		PVT->GetBw(n_act_cell_, pw_, bw_);	//Calculate Bw at time step n+1
		PVT->GetVisO(n_act_cell_, po_, vo_); //Calculate Vo at time step n+1
		PVT->GetVisW(n_act_cell_, pw_, vw_); //Calculate Vw at time step n+1
		PVT->GetDenO(n_act_cell_, bo_, deno_); //Calculate Density of Oil at time step n+1
		PVT->GetDenW(n_act_cell_, bw_, denw_); //Calculate Density of Water at time step n+1
		SAT->GetKro(n_act_cell_, so_, kro_); //Calculate Kro at time step n+1
		SAT->GetKrw(n_act_cell_, sw_, krw_); //Calculate Krw at time step n+1

		// Assemble mobility, accumulation term
		CalMobiOil();  //Calculate mobility of oil at time step n+1
		CalMobiWater(); //Calculate mobility of water at time step n+1

		CalAccuOil(Grid, dt); //Calculate accumulation of oil  
		CalAccuWater(Grid, dt); //Calculate accumulation of water 

		// Calculate Connection-wise propeties
		CalGammaOil(); //Calculate gamma of oil at all connections
		CalGammaWater(); //Calculate gamma of water at all connections
		CalDPhiOil(Grid); //Calculate potential difference of oil at all connections
		CalDPhiWater(Grid); //Calculate potential difference of water at all connections

		//Assemble flux term
		CalFo(); //Calculate flux term of oil at all connections
		CalFw(); //Calculate flux term of water at all connections

	}
	else if (phase_type_ == SimCtrl::OWG){

	}
}

double CState::CalResidualNormInf(CGrid *Grid, double dt){
	//	cout << "dt:" << dt << endl;
	double norm_resid = -1;
	double normalized_resid_w, normalized_resid_o;
	const double *v = Grid->GetVolume();
	//	cout<<"V:"<< v[0] <<" Poro:"<<poro[0]<<" bw:"<<bw_[0]<<endl;
	for (int n = 0; n < num_phase_*n_act_cell_; n++){
		if (n % 1 == 0){
			// n is odd, correspond to water
			normalized_resid_w = 5.615 / (v[n / 2] * poro_[n / 2])* bw_[n / 2 ]* dt*resid_[n];
			if (fabs(normalized_resid_w)>norm_resid){ 
				norm_resid = fabs(normalized_resid_w); 
				//				cout <<"Norm: "<< norm_resid <<" "<<fabs(normalized_resid_w) << endl;
			}
		}
		else{
			normalized_resid_o = 5.615 / (v[(n - 1) / 2] * poro_[(n - 1) / 2])* bo_[(n - 1) / 2] * dt*resid_[n];
			if (fabs(normalized_resid_o)>norm_resid){ norm_resid = fabs(normalized_resid_o); }
		}
	}
	//	cout << norm_resid << endl;
	return norm_resid;
}

double CState::AssembleResidual(CGrid* Grid, double dt, vector<CStandardWell*> &std_well){
	double err;
	switch (phase_type_)
	{
	case SimCtrl::OWG:
		// To be completed
		break;
	case SimCtrl::OW:
		//		cout << "Assemble residual OW..." << endl;
		//Assemble Residual Vector
		AssembleResidualOW(std_well);
		//		Residual_Debug("residual_after_assemble.txt");
		err = CalResidualNormInf(Grid,dt);
		break;
	case SimCtrl::OG:
		//===========CalDPoroDPo==============================================
		//Assemble Residual Vector
		AssembleResidualOG(std_well);
		err = CalResidualNormInf(Grid, dt);
		break;
	default:
		break;
	}
	return err;
}

bool CState::AssembleJacobian(CGrid* Grid, double dt, vector<CStandardWell*> &std_well){
	switch (phase_type_)
	{
	case SimCtrl::OWG:
		// To be completed
		break;
	case SimCtrl::OW:
		//Assemble Jacobian Matrix
		AssembleJacobianOW(Grid, dt, std_well);
		break;
	case SimCtrl::OG:
		//Assemble Jacobian Matrix
		AssembleJacobianOG(Grid, dt, std_well);
		break;
	default:
		break;
	}
	return true;
}

void CState::OutPutFPR(CSchedule &SCH){
	double fpr=0;
	for (int i = 0; i < n_act_cell_; i++){
		fpr += po_[i];
	}
	fpr = fpr / n_act_cell_;
	fstream fout("FPR.txt", ios::out | ios::app);
	fout << SCH.GetTCurrent() << "  " << fpr << endl;
}

void CState::InputData(){
	string fn1("p_n.dat");
	ifstream in;
	in.open(fn1.c_str());
	for (int i = 0; i < n_act_cell_; i++){
		in >> po_n_[i];
		pw_n_[i] = po_n_[i];
	}
	in.close();

	string fn2("sw_n.dat");
	in.open(fn2.c_str());
	for (int i = 0; i < n_act_cell_; i++){
		in >> sw_n_[i];
		so_n_[i] = 1 - sw_n_[i];
	}
	in.close();

	string fn3("p.dat");
	in.open(fn3.c_str());
	for (int i = 0; i < n_act_cell_; i++){
		in >> po_[i];
		pw_[i] = po_[i];
	}
	in.close();

	string fn4("sw.dat");
	in.open(fn4.c_str());
	for (int i = 0; i < n_act_cell_; i++){
		in >> sw_[i];
		so_[i] = 1 - sw_[i];
	}
	in.close();
}

void CState::OutputJacobian(){
	//
	string fn1("jaco_d.dat");
	ofstream out;
	out.open(fn1.c_str());
	out.precision(15);
	for (int i = 0; i < n_act_cell_; i++){
		out << jaco_d_[i](0, 0) << "  " << jaco_d_[i](0, 1) << endl;
		out << jaco_d_[i](1, 0) << "  " << jaco_d_[i](1, 1) << endl;
		out << endl;
	}
	out.close();
	//
	string fn2("jaco_ud.dat");
	out.open(fn2.c_str());
	out.precision(15);
	for (int i = 0; i < nconns_; i++){
		out << jaco_ud_[i](0, 0) << "  " << jaco_ud_[i](0, 1) << endl;
		out << jaco_ud_[i](1, 0) << "  " << jaco_ud_[i](1, 1) << endl;
		out << endl;
	}
	out.close();
	//
	string fn3("jaco_ld.dat");
	out.open(fn3.c_str());
	out.precision(15);
	for (int i = 0; i < nconns_; i++){
		out << jaco_ld_[i](0, 0) << "  " << jaco_ld_[i](0, 1) << endl;
		out << jaco_ld_[i](1, 0) << "  " << jaco_ld_[i](1, 1) << endl;
		out << endl;
	}
	out.close();
	// Output Residual
	string fn4("residual.dat");
	out.open(fn4.c_str());
	out.precision(15);
	for (int i = 0; i < num_phase_*n_act_cell_; i++){
		out <<resid_[i] << endl;
	}
	out.close();
}

double CState::CalCFLOW(CGrid* Grid,double dt){
	double cfl = -1.0;
	double temp = 0.0;
	const double* v = Grid->GetVolume();
	double* outflux = new double[n_act_cell_];
	memset(outflux, 0, n_act_cell_*sizeof(double));
	int l, r;
	for (int i = 0; i < nconns_; i++){
		l = l_[i];
		r = r_[i];
		if (fo_[i]>0){ //Flows from l to r
			outflux[l] += fo_[i];
		}
		else{
			outflux[r] -= fo_[i];
		}
		if (fw_[i] > 0){
			outflux[l] += fw_[i];
		}
		else{
			outflux[r] -= fw_[i];
		}
	}
	for (int i = 0; i < n_act_cell_; i++){
		temp = 5.615* outflux[i] * dt / (poro_[i] * v[i]);
		if (temp>cfl){
			cfl = temp;
		}
	}
	delete[] outflux;
	return cfl;
}
