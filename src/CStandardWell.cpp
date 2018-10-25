#include "CStandardWell.h"
#include "CGrid.h"
#include <math.h>

//=====================================================
CStandardWell::CStandardWell(string well_name, int n){
	well_name_ = well_name;
	n_ = n;
	//Set up default values
	po_ = 0.0;
	pg_ = 0.0;
	pw_ = 0.0;
	mo_ = 0.0;
	mg_ = 0.0;
	mw_ = 0.0;
	bo_ = 0.0;
	bg_ = 0.0;
	bw_ = 0.0;
	rso_ = 0.0;
	dmo_dp_ = 0.0;
	dmg_dp_ = 0.0;
	dmw_dp_ = 0.0;
	dmo_dsw_ = 0.0;
	dmw_dsw_ = 0.0;
	dmg_dsw_ = 0.0;
	dmo_dsg_ = 0.0;
	dmg_dsg_ = 0.0;
	dqo_dp_ = 0.0;
	dqo_dsw_ = 0.0;
	dqo_dsg_ = 0.0;
	dqw_dp_ = 0.0;
	dqw_dsw_ = 0.0;
	dqw_dsg_ =0.0;
	dqg_dp_ = 0.0;
	dqg_dsw_ = 0.0;
	dqg_dsg_ = 0.0;
	//
	ctrl_mode_ = CBHP;
	TL_ORAT_ = 1e20;
	TL_BHP_ = 0.0;
	TL_WRAT_ = 1e20;
	TL_GRAT_ = 1e20;
	TL_LRAT_ = 1e20;
	//
	ORAT_ = 0.0;
//	ORAT_RC = 0.0;
	GRAT_ = 0.0;
	WRAT_ = 0.0;
//	WRAT_RC = 0.0;
	BHP_ = 0.0;
	LRAT_ = 0.0;

}
CStandardWell::~CStandardWell(){
	
}
bool CStandardWell::set_r(double r){
	r_ = r;
	return true;
}
string CStandardWell::get_well_name() const{
	return well_name_;
}
bool CStandardWell::CalWellIndex(CGrid* grid){
	double kx, ky, dx, dy,dz;
	kx = grid->GetKx(n_);
	ky = grid->GetKy(n_);
	dx = grid->GetDX(n_);
	dy = grid->GetDY(n_);
	dz = grid->GetDZ(n_);
	r0_ = 0.28  *(pow((pow(ky / kx, 0.5)*dx*dx + pow(kx / ky, 0.5)*dy*dy), 0.5)) / (pow(ky / kx, 0.25) + pow(kx / ky, 0.25));
	WI_ = alpha * 2 * PI*pow(kx*ky, 0.5)*dz / log( r0_ / r_);
	return true;
}
int CStandardWell::get_n(){
	return n_;
}
bool CStandardWell::set_well_index(double WI){
	WI_ = WI;
	return true;
}
//=====================================================
double CStandardWell::get_qo(){
	return ORAT_;
}

double CStandardWell::get_qw(){
	return WRAT_;
}

double CStandardWell::get_qg(){
	return GRAT_;
}

double CStandardWell::get_dqo_dp(){
	return dqo_dp_;
}
double CStandardWell::get_dqw_dp(){
	return dqw_dp_;
}
double CStandardWell::get_dqg_dp(){
	return dqg_dp_;
}
double CStandardWell::get_dqo_dsw(){
	return dqo_dsw_;
}
double CStandardWell::get_dqw_dsw(){
	return dqw_dsw_;
}
double CStandardWell::get_dqg_dsw(){
	return dqg_dsw_;
}
double CStandardWell::get_dqo_dsg(){
	return dqo_dsg_;
}
double CStandardWell::get_dqw_dsg(){
	return dqw_dsg_;
}
double CStandardWell::get_dqg_dsg(){
	return dqg_dsg_;
}

//=====================================================
bool  CStandardWell::set_ctrl_mode(CTRLMODE ctrl_mode){
	ctrl_mode_ = ctrl_mode;
	return true;
}

bool CStandardWell::set_TL_ORAT(double TL_ORAT){
	TL_ORAT_ = TL_ORAT;
	return true;
}

bool CStandardWell::set_TL_BHP(double TL_BHP){
	TL_BHP_ = TL_BHP;
	return true;
}

bool CStandardWell::set_TL_WRAT(double TL_WRAT){
	TL_WRAT_ = TL_WRAT;
	return true;
}

bool CStandardWell::set_TL_LRAT(double TL_LRAT){
	TL_LRAT_ = TL_LRAT;
	return true;
}

bool CStandardWell::set_TL_GRAT(double TL_GRAT){
	TL_GRAT_ = TL_GRAT;
	return true;
}
//=================================================
bool CStandardWell::insert_multistage_control(double start_time, CTRLMODE ctrl_mode, double target){
	if(multistage_start_time.empty()){
		multistage_start_time.push_back(start_time);
		multistage_control_mode.push_back(ctrl_mode);
		multistage_target.push_back(target);
	}else{
		for(int k =0; k < multistage_start_time.size(); k++){
			if (start_time < multistage_start_time[k]){
				multistage_start_time.insert(multistage_start_time.begin()+k, start_time);
				multistage_control_mode.insert(multistage_control_mode.begin()+k, ctrl_mode);
				multistage_target.insert(multistage_target.begin()+k, target);
				break;
			}else if(start_time == multistage_start_time[k]){
				break;
			}
			if(k == (multistage_start_time.size()-1)){
				multistage_start_time.push_back(start_time);
				multistage_control_mode.push_back(ctrl_mode);
				multistage_target.push_back(target);
			}
		}
	}
	return true;
}


bool CStandardWell::update_control(double t){
	double eps = 0.00001;
	string control_mode;
	if(!multistage_start_time.empty()){
		if(t >= multistage_start_time[0] || fabs(t - multistage_start_time[0]) < eps){

			set_ctrl_mode(multistage_control_mode[0]);

			switch(multistage_control_mode[0])	{
			case CStandardWell::CBHP:
				set_TL_BHP(multistage_target[0]);
				BHP_ = multistage_target[0];
				control_mode = "BHP";
				break;
			case CStandardWell::CORAT:
				set_TL_ORAT(multistage_target[0]);
				ORAT_ = multistage_target[0];
				control_mode = "ORAT";
				break;
			case CStandardWell::CWRAT:
				set_TL_WRAT(multistage_target[0]);
				WRAT_ = multistage_target[0];
				control_mode = "WRAT";
				break;
			case CStandardWell::CLRAT:
				set_TL_LRAT(multistage_target[0]);
				LRAT_ = multistage_target[0];
				control_mode = "LRAT";
				break;
			case CStandardWell::CGRAT:
				set_TL_GRAT(multistage_target[0]);
				GRAT_ = multistage_target[0];
				control_mode = "GRAT";
				break;
			default:
				break;
			}
			cout<< well_name_ << " switch to constant "  << control_mode << " of " << multistage_target[0] << " at day " << t << endl;
			return true;
		}
		return false;
	}
	return false;
}

bool CStandardWell::remove_control_stage(double t){
	double eps = 0.00001;
	if(!multistage_start_time.empty()){
		if (t >= multistage_start_time[0] || fabs(t - multistage_start_time[0]) < eps){
//		if (t >= multistage_start_time[0]){
			multistage_start_time.erase(multistage_start_time.begin());
			multistage_control_mode.erase(multistage_control_mode.begin());
			multistage_target.erase(multistage_target.begin());
		}
	}
	return true;
}

//=================================================
bool  CStandardWell::RecordResult(double t){
	TIME_.push_back(t);
	ORAT_REC_.push_back(ORAT_);
	WRAT_REC_.push_back(WRAT_);
	GRAT_REC_.push_back(GRAT_);
	LRAT_REC_.push_back(ORAT_ + WRAT_ + GRAT_);
	BHP_REC_.push_back(BHP_);
	WWCT_REC_.push_back(ORAT_/(ORAT_+WRAT_+GRAT_));
	WBP_REC_.push_back(po_);
	return true;
}

vector<double> CStandardWell::get_TIME() const{
	return TIME_;
}
vector<double> CStandardWell::get_ORAT() const{
	return ORAT_REC_;
}
vector<double> CStandardWell::get_WRAT() const{
	return WRAT_REC_;
}
vector<double> CStandardWell::get_BHP() const{
	return BHP_REC_;
}
vector<double> CStandardWell::get_LRAT() const{
	return LRAT_REC_;
}
vector<double> CStandardWell::get_GRAT() const{
	return GRAT_REC_;
}

vector<double> CStandardWell::get_WBP() const{
	return WBP_REC_;
}

vector<double> CStandardWell::get_WWCT() const{
	return WWCT_REC_;
}
