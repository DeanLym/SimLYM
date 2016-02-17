/*
 * CHistoryMatching.hpp
 *
 *  Created on: Feb 16, 2016
 *      Author: yiminliu
 */

#ifndef CHISTORYMATCHING_HPP_
#define CHISTORYMATCHING_HPP_

#include<vector>
#include<string>

using namespace std;

class CHMTarget{
public:
	enum DATA_TYPE {ORAT,WRAT,WCT,BHP};
	enum DATA_DOMAIN {FIELD,GROUP,WELL};
	enum STD_TYPE {FRACTION,ABSOLUTE};
	CHMTarget(DATA_DOMAIN domain,string wd_name,DATA_TYPE data_type,double t, double value, double std);
public:
	CHMTarget& operator=(const CHMTarget &other);
public:
	double t_;
	double value_;
	double std_;
	DATA_DOMAIN domain_;
	string wd_name_;
	DATA_TYPE data_type_;
	STD_TYPE std_type_;
	double std_min_;
	double std_max_;
};

class CStandardWell;

class CHistoryMatching{
public:
	void SetHMTarget(vector<CHMTarget> hm_target);
	void SetHMTarget(char * hist_file);
	void cal_d_sim(const vector<CStandardWell *> well_list);
	void InputHist(vector<double> d_obs, vector<double> std_obs);
	void OutputHist(vector<double> d_obs, vector<double> std_obs);
	double GetDataMismatch(const vector<CStandardWell *> well_list);
	CHMTarget::DATA_TYPE Str2DataType(string str);
	CHMTarget::DATA_DOMAIN Str2DataDomain(string str);
protected:
	double GetWellData(const vector<CStandardWell *> well_list, string name,CHMTarget::DATA_TYPE data_type,double t);
private:
	int num_data_;
	vector<CHMTarget> hm_target_;
	vector<double> d_sim_;
};






#endif /* CHISTORYMATCHING_HPP_ */
