#ifndef __CGrid_H__
#define __CGrid_H__

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "CConn.h"
#include <fstream>
#include <string.h>
#include <iostream>

using namespace std;

class CGrid{
public:
	//Constants
	static constexpr double beta = 0.006944;
	static constexpr double alpha = 0.001127;

public:

	CGrid(int ncell){
		ncell_ = ncell;
		poro_= new double[ncell_];
		memset(poro_, 0, ncell_*sizeof(double));
		poro_n_= new double[ncell_];
		memset(poro_n_, 0, ncell_*sizeof(double));
		v_ = new double[ncell_];
		memset(v_, 0, ncell_*sizeof(double));
		d_ = new double[ncell_];
		memset(d_, 0, ncell_*sizeof(double));
		dmin_ = 1e8;
		dmax_ = -1e8;
		cr_ = 0;
		p_ref_ = 1;
	}

	int get_ncell(){
		return ncell_;
	}

	int get_nconns(){
		return nconns_;
	}

	bool InputPoro(double poro){
		for (int i = 0; i < ncell_; i++){
			poro_[i] = poro;
		}
		return true;
	}

	bool InputCr(double cr){
		cr_ = cr;
		return true;
	}

	bool InputPref(double p_ref){
		p_ref_ = p_ref;
		return true;
	}

	bool InputPoro(double* poro){
		for (int i = 0; i < ncell_; i++){
			poro_[i] = poro[i];
		}
		return true;
	}

	bool InputVolume(double v){
		for (int i = 0; i < ncell_; i++){
			v_[i] = v;
		}
		return true;
	}

	bool InputVolume(double* v){
		for (int i = 0; i < ncell_; i++){
			v_[i] = v[i];
		}
		return true;
	}
	bool InputVolume(char* file){
		ifstream in;
		in.open(file);
		if (in.is_open()){
			for (int i = 0; i < ncell_; i++){
				in >> v_[i];
			}
		}
		else{
			cout << "Open Volume File Failed!" << endl;
			exit(1);
			return false;
		}
		in.close();
		return true;
	}

	bool InputDepth(double* d){
		for (int i = 0; i < ncell_; i++){
			d_[i] = d[i];
		}
		return true;
	}

	bool InputDepth(double d){
		for (int i = 0; i < ncell_; i++){
			d_[i] = d;
		}
		return true;
	}

	double GetPoro(int k){
		return poro_[k];
	}

    const double* GetPoro(){
    	return poro_;
    }

    bool GetPoro(double *po, double* poro){
    	double dp;
    	for(int i=0; i < ncell_; i++){
    		dp = po[i] - p_ref_;
    		poro[i] = poro_[i]*(1+cr_*dp+0.5*cr_*cr_*dp*dp);
    	}
    	return true;
    }

    bool GetDPoroDPo(double *po, double *dporo_dpo){
    	for(int i=0; i < ncell_; i++){
    		dporo_dpo[i] = poro_[i]*cr_*(1+cr_*(po[i] - p_ref_));
    	}
    	return true;
    }

	double GetVolume(int k){
		return v_[k];
	}

    const double* GetVolume(){
        return v_;
    }

	double GetDepth(int k){
		return d_[k];
	}

    const double* GetDepth(){
        return d_;
    }

	bool InputConnList(ConnList* conn_list,int nconns){
		conn_list_ = conn_list;
		for (int k = 0; k < nconns; k++){
			conn_list_[k] = conn_list[k];
		}
		return true;
	}

	bool InputConnList ( int *l, int *r, double *trans, int nconns){
		nconns_ = nconns;
        conn_list_ = new ConnList(nconns_);
        conn_list_->InsertConnection(l,r,trans);
		return true;
	}

	bool InputConnList(int nconns, char* file){
		nconns_ = nconns;
		conn_list_ = new ConnList(nconns_);

		ifstream in;
		in.open(file);
		if (in.is_open()){
			int *l = new int[nconns_];
			int *r = new int[nconns_];
			double *trans = new double[nconns_];
			for (int i = 0; i < nconns_; i++){
				in >> l[i];
				in >> r[i];
				in >> trans[i];
			}
			conn_list_->InsertConnection(l, r, trans);
			delete[] l;
			delete[] r;
			delete[] trans;
		}
		else{
			cout << "Open Connectivity List Failed!" << endl;
			exit(1);
			return false;
		}
		in.close();
		return true;
	}

    const ConnList* GetConnList(){
        return conn_list_;
    }
	bool CalDeltaD(){
		delta_d_ = new double[nconns_];
		memset(delta_d_, 0, nconns_*sizeof(double));
        const int* l = conn_list_->get_l();
        const int* r = conn_list_->get_r();
		for (int k = 0; k < nconns_; k++){
            delta_d_[k] = d_[l[k]] - d_[r[k]];
		}
		return true;
	}
	double GetdMax(){
		return dmax_;
    }

	double GetdMin(){
		return dmin_;
	}

    const double* GetDeltaD(){
        return delta_d_;
    }
	bool OutputConnList(){
		string fn("connlist.dat");
		ofstream out;
		out.open(fn.c_str());
		const int *l = conn_list_->get_l();
		const int *r = conn_list_->get_r();
		const double* trans = conn_list_->get_trans();
		for (int i = 0; i < nconns_; i++){
			out << l[i] << "  " << r[i] << "  " << trans[i] << endl;
		}
		out.close();
		return true;
	}
public:
	// Virtual functions for derived classes
	// For structured grid
	virtual bool InputKx(double kx) = 0;
	virtual bool InputKx(double *kx) = 0;
	virtual bool InputKx(char *file_name) = 0;
	virtual bool InputKy(double ky) = 0;
	virtual bool InputKy(double *ky) = 0;
	virtual bool InputKy(char *file_name) = 0;
	virtual bool InputKz(double kz) = 0;
	virtual bool InputKz(double *kz) = 0;
	virtual bool InputKz(char *file_name) = 0;
	virtual int GetnX() = 0;
	virtual int GetnY() = 0;
	virtual int GetnZ() = 0;
	virtual bool InputDx(double dx) = 0;
	virtual bool InputDx(double *dx) = 0;
	virtual bool InputDy(double dy) = 0;
	virtual bool InputDy(double *dy) = 0;
	virtual bool InputDz(double dz) = 0;
	virtual bool InputDz(double *dz) = 0;
	virtual double GetDX(int n) = 0;
	virtual double GetDY(int n) = 0;
	virtual double GetDZ(int n) = 0;
	virtual double GetKx(int n) = 0;
	virtual double GetKy(int n) = 0;
	virtual double GetKz(int n) = 0;
	// For Cartesian Grid
	virtual bool CalVolume() = 0;
	virtual bool CalDepth() = 0;
	virtual bool GenerateConnList() = 0;

	virtual int GetIndex(int i, int j, int k) = 0;
	virtual bool InputTops(double tops) = 0;
public:
	~CGrid(){
		delete[] poro_;
		delete[] poro_n_;
		delete[] v_;
		delete[] delta_d_;
		delete[] d_;
		delete conn_list_;
	}
protected:
	int ncell_;
	int nconns_;
//	vector<double> poro_;
	double* poro_;
	double* poro_n_;
	double p_ref_;
	double cr_;

	//	vector<double> v_;
	double* v_;
//	vector<double> d_;
	double* d_;
	double* delta_d_;

	ConnList* conn_list_;
	double dmin_;
	double dmax_;


};

#endif
