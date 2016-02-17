
#include "CHistoryMatching.hpp"
#include "CStandardWell.h"
#include <exception>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <fstream>

using namespace std;
//=====================CHistoryMatching========================

void CHistoryMatching::SetHMTarget(vector<CHMTarget> hm_target){
	hm_target_ = hm_target;
	num_data_ = hm_target_.size();
}

CHMTarget::DATA_TYPE CHistoryMatching::Str2DataType(string str){
	if(str == "ORAT")
		return CHMTarget::ORAT;
	if(str == "WRAT")
		return CHMTarget::WRAT;
	if(str == "WCT")
		return CHMTarget::WCT;
	if(str == "BHP")
		return CHMTarget::BHP;
	throw runtime_error("Unknown data type string ...");
	return  CHMTarget::ORAT;
}

CHMTarget::DATA_DOMAIN  CHistoryMatching::Str2DataDomain(string str){
	if(str == "WELL")
		return CHMTarget::WELL;
	if(str == "GROUP")
		return CHMTarget::GROUP;
	if(str == "FIELD")
		return CHMTarget::FIELD;
	throw runtime_error("Unknown data domain string ...");
	return CHMTarget::WELL;
}

void CHistoryMatching::SetHMTarget(char* hist_file){
	ifstream in(hist_file);
	if(!in.is_open())
		throw runtime_error("Can not open hist file...");
	in >> num_data_;
	CHMTarget::DATA_DOMAIN data_domain;
	CHMTarget::DATA_TYPE data_type;
	string wd_name;
	string temp1,temp2;
	double t, d_obs, std_obs;
	hm_target_.clear();
	for(int i=0;i<num_data_;i++){
		in >> temp1 >> wd_name >> temp2 >> t >> d_obs >> std_obs;
		data_domain = Str2DataDomain(temp1);
		data_type = Str2DataType(temp2);
		hm_target_.push_back(CHMTarget(data_domain,wd_name,data_type,t,d_obs,std_obs));
	}
	in.close();
}

double CHistoryMatching::GetDataMismatch(const vector<CStandardWell *> well_list){
	cal_d_sim(well_list);
	double Sd = 0.0, temp;
	for(int i=0;i<hm_target_.size();i++){
		temp = (d_sim_[i] - hm_target_[i].value_)/hm_target_[i].std_;
		Sd += temp * temp;
	}
	return Sd;
}


void CHistoryMatching::cal_d_sim(const vector<CStandardWell *> well_list){
	for(int i=0;i<num_data_;i++){
		switch(hm_target_[i].domain_){
		case CHMTarget::WELL:
			d_sim_.push_back(
					GetWellData(well_list,
							hm_target_[i].wd_name_,
							hm_target_[i].data_type_,
							hm_target_[i].t_));
			break;
		case CHMTarget::GROUP:
			// To be complete
		case CHMTarget::FIELD:
			// To be complete
		default:
			throw runtime_error("Unknown domain of history matching target...");
		}
	}
	//	for(int i=0;i<num_data_;i++){
	//		double std;
	//		switch(hm_target_[i].std_type_){
	//		case CHMTarget::FRACTION:
	//			std = d_obs_[i] * hm_target_[i].std_;
	//			break;
	//		case CHMTarget::ABSOLUTE:
	//			std = hm_target_[i].std_;
	//			break;
	//		default:
	//			throw runtime_error("Unknown STD Type...");
	//			break;
	//		}
	//		std = std>hm_target_[i].std_min_?std:hm_target_[i].std_min_;
	//		std = std<hm_target_[i].std_max_?std:hm_target_[i].std_max_;
	//		std_obs_.push_back(std);
	//	}
}


//void CHistoryMatching::InputHist(vector<double> d_obs, vector<double> std_obs){
//	d_obs_ = d_obs;
//	std_obs_ = std_obs;
//}

struct dbl_cmp {
	dbl_cmp(double v, double d):val(v), delta(d) { }
	inline bool operator()(const double &x) const {
		return fabs(x-val) < delta;
	}
private:
	double val, delta;
};


double CHistoryMatching::GetWellData(const vector<CStandardWell *> well_list, string name,CHMTarget::DATA_TYPE data_type,double t){
	vector<double> data_all_time;
	vector<double> T;
	int name_in_list = 0;
	for(int i=0;i<well_list.size();i++){
		if(well_list[i]->get_well_name() == name){
			name_in_list = 1;
			T = well_list[i]->get_TIME();
			switch(data_type){
			case CHMTarget::ORAT:
				data_all_time = well_list[i]->get_ORAT();
				break;
			case CHMTarget::WRAT:
				data_all_time = well_list[i]->get_WRAT();
				break;
			case CHMTarget::WCT:
				data_all_time = well_list[i]->get_WWCT();
				break;
			case CHMTarget::BHP:
				data_all_time = well_list[i]->get_BHP();
				break;
			default:
				throw runtime_error("Unknown data type of history matching target...");
			}
		}
	}
	if(name_in_list == 0)
		throw runtime_error("Can not find well name: " + name);

	double epi = 1e-3;
	auto ind = std::find_if(T.begin(),T.end(),dbl_cmp(t,epi));
	if(ind != T.end())
		return data_all_time[ind - T.begin()];
	else
		throw runtime_error("Can not find data at specified time..");
	return 0.0;
}



//=====================CHMTarget========================


CHMTarget::CHMTarget(DATA_DOMAIN domain,string wd_name,DATA_TYPE data_type,double t,double value, double std):
		domain_(domain),
		wd_name_(wd_name),
		data_type_(data_type),
		t_(t),
		value_(value),
		std_type_(FRACTION),
		std_min_(-1.0),
		std_max_(1e20),
		std_(std)
{

}


CHMTarget& CHMTarget::operator=(const CHMTarget &other){
	this->data_type_ = other.data_type_;
	this->domain_ = other.domain_;
	this->wd_name_ = other.wd_name_;
	this->data_type_ = other.data_type_;
	return *this;
}











