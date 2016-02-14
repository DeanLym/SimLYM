#ifndef __CSAT_H__
#define __CSAT_H__

#include <stdio.h>
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

class CSAT{
public:
	CSAT():
		sw_tab_(NULL),
		so_tab_(NULL),
		sg_tab_(NULL),
		krw_tab_(NULL),
		krg_tab_(NULL),
		krog_tab_(NULL),
		krow_tab_(NULL),
		dkrw_dsw_tab_(NULL),
		dkrow_dsw_tab_(NULL),
		dkrg_dsg_tab_(NULL),
		dkrog_dsg_tab_(NULL)
		{

		}
	~CSAT(){
		delete[] sw_tab_; sw_tab_ = NULL;
		delete[] krw_tab_; krw_tab_ = NULL;
		delete[] so_tab_; so_tab_ = NULL;
		delete[] krow_tab_; krow_tab_ = NULL;
		delete[] dkrow_dsw_tab_; dkrow_dsw_tab_ = NULL;
		delete[] dkrw_dsw_tab_; dkrw_dsw_tab_ = NULL;
		delete[] sg_tab_; sg_tab_ = NULL;
		delete[] krg_tab_; krg_tab_ = NULL;
		delete[] krog_tab_; krog_tab_ = NULL;
		delete[] dkrg_dsg_tab_; dkrg_dsg_tab_ = NULL;
		delete[] dkrog_dsg_tab_; dkrog_dsg_tab_ = NULL;
	}

	bool GetKrw(int size,double *sw,double *krw){
		for(int i=0;i<size;i++){
			krw[i] = PiecewiseLinearInter(nr_sat_tab_, sw_tab_, krw_tab_, sw[i]);
		}

		return true;
	}
	bool GetKro(int size,double *so,double *kro){
		for(int i=0;i<size;i++){
			kro[i] = PiecewiseLinearInter(nr_sat_tab_, sw_tab_, krow_tab_, 1 - so[i]);
		}

		return true;
	}

	double GetKrwPDSw(int size,double *sw,double *dkrw_dsw){
		for(int i=0;i<size;i++)
			dkrw_dsw[i] = PiecewiseLinearInter(nr_sat_tab_, sw_tab_, dkrw_dsw_tab_, sw[i]);
		return true;
	}

	bool GetKroPDSw(int size,double *sw, double *dkrow_dsw){
		for(int i=0;i<size;i++)
			dkrow_dsw[i] = PiecewiseLinearInter(nr_sat_tab_, sw_tab_, dkrow_dsw_tab_, sw[i]);
		return true;
	}

	void SetSoEndPoints(double sorg){
		sorg_ = sorg;
	}
	void SetSgEndPoints(double sgco, double sgcr){
		sgco_ = sgco;
		sgcr_ = sgcr;
	}
	double Getsgco_(){ return sgco_; }
	double Getsgcr_(){ return sgcr_; }
	void SetSwEndPoints(double swco){
		swco_ = swco;
	}
	double Getswco_(){ return swco_; }
	double GetKroPDSg(double sg){
		double res=0.0;
		res = PiecewiseLinearInter(nr_sat_tab_, sg_tab_, dkrog_dsg_tab_, sg);
		return res;
	}

	bool GetKroPDSg(int size, double *sg, double *dkrog_dsg){
		for(int i=0;i<size;i++)
			dkrog_dsg[i] = PiecewiseLinearInter(nr_sat_tab_, sg_tab_, dkrog_dsg_tab_, sg[i]);
		return true;
	}
	bool GetKrgPDSg(int size,double *sg, double *dkrg_dsg){
		for(int i=0;i<size;i++)
			dkrg_dsg[i] = PiecewiseLinearInter(nr_sat_tab_, sg_tab_, dkrg_dsg_tab_, sg[i]);
		return true;
	}
	bool GetKrgPDSw(int size,double *sw,double* dkrg_dsw){
		for(int i=0;i<size;i++)
			dkrg_dsw[i] = 0.0;
		return true;
	}
	bool GetKrg(int size, double *sg, double *krg){
		for (int i = 0; i<size; i++)
			krg[i]=PiecewiseLinearInter(nr_sat_tab_, sg_tab_, krg_tab_, sg[i]);
		return true;
	}
	void SetSWOF(int _nr_sat_tab_, char *satFile){
		nr_sat_tab_ = _nr_sat_tab_;
		sw_tab_ = new double[nr_sat_tab_];
		krw_tab_ = new double[nr_sat_tab_];
		so_tab_ = new double[nr_sat_tab_];
		krow_tab_ = new double[nr_sat_tab_];
		dkrow_dsw_tab_ = new double[nr_sat_tab_];
		dkrw_dsw_tab_ = new double[nr_sat_tab_];

		fstream fp;
		fp.open(satFile, ios::in);

		string temp;
		fp >> temp;
		double dummy_capillary;
		for (int i = 0; i < nr_sat_tab_; i++){
			fp >> sw_tab_[i];
			so_tab_[i] = 1 - sw_tab_[i];
			fp >> krw_tab_[i];
			cout << krw_tab_[i] << endl;
			fp >> krow_tab_[i];
			fp >> dummy_capillary;
		}
		fp.close();


		CreatePDTable(nr_sat_tab_,krw_tab_,sw_tab_,dkrw_dsw_tab_);
		CreatePDTable(nr_sat_tab_,krow_tab_,sw_tab_,dkrow_dsw_tab_);
	}
	void SetSGOF(int _nr_sat_tab_,char *satFile){
		nr_sat_tab_ = _nr_sat_tab_;
		sg_tab_ = new double[nr_sat_tab_];
		krg_tab_ = new double[nr_sat_tab_];
		so_tab_ = new double[nr_sat_tab_];
		krog_tab_ = new double[nr_sat_tab_];
		dkrog_dsg_tab_ = new double[nr_sat_tab_];
		dkrg_dsg_tab_ = new double[nr_sat_tab_];


		fstream fp;
		fp.open(satFile, ios::in);

		string temp;
		fp >> temp;
		double dummy_capillary;
		for (int i = 0; i < nr_sat_tab_; i++){
			fp >> sg_tab_[i];
			so_tab_[nr_sat_tab_-i-1] = 1 - sg_tab_[i];
			fp >> krg_tab_[i];
			fp >> krog_tab_[nr_sat_tab_ - i - 1];
			fp >> dummy_capillary;
		}
		fp.close();

		CreatePDTable(nr_sat_tab_,krg_tab_,sg_tab_,dkrg_dsg_tab_);
		CreatePDTable(nr_sat_tab_,krog_tab_,sg_tab_,dkrog_dsg_tab_);

	}
	double PiecewiseLinearInter(int n, double *pX, double *pY, double x){
		double ans = 0.0;
		if(pX[0]>=x){
			ans = pY[0];
			return ans;
		}else if(pX[n-1]<=x){
			ans = pY[n-1];
			return ans;
		}
		else if (pX[0] <= x&&pX[1] > x){
			ans = (pY[1] - pY[0]) / (pX[1] - pX[0])*(x - pX[0]) + pY[0];
			return ans;
		}
		else if (pX[n-2] <= x&&pX[n-1] >= x){
			ans = (pY[n-1] - pY[n-2]) / (pX[n-1] - pX[n-2])*(x - pX[n-2]) + pY[n-2];
			return ans;
		}
		else{
			for (int i = 1; i < n - 1; i++){
				if (pX[i] <= x&&pX[i + 1]>x){
					ans = (pY[i + 1] - pY[i]) / (pX[i + 1] - pX[i])*(x - pX[i]) + pY[i];
					return ans;
				}
			}
		}
	}
	bool CreatePDTable(int n,double *tab_y, double* tab_x, double *tab_dy_dx){
		tab_dy_dx[0] = (tab_y[1] - tab_y[0])/(tab_x[1]-tab_x[0]);
		for(int i=1 ; i<n -1 ; i++){
			tab_dy_dx[i] = (tab_y[i+1] - tab_y[i-1])/(tab_x[i+1]-tab_x[i-1]);
		}
		tab_dy_dx[n-1] =
				(tab_y[n - 1 ] - tab_y[n - 2])
				/(tab_x[n - 1 ]-tab_x[n - 2]);
		return true;
	}
	//	double GetPcOW(double sw);
	//	double GetPcOG(double sg);
private:
	double swcr_, swco_, sw_max_;
	double sorg_;
	double sgcr_, sgco_;
	double *sg_tab_, *krg_tab_;
	double *so_tab_, *krog_tab_;
	double *sw_tab_, *krw_tab_, *krow_tab_;
	double *dkrow_dsw_tab_, *dkrw_dsw_tab_;
	double *dkrog_dsg_tab_, *dkrg_dsg_tab_;
	int nr_sat_tab_;
	//	double *pSw, *pSo, *pSg;
	//	double *pKrw, *Kro, *pKrg;
	//	double *pPcOW, *pPcOG;
};

#endif
