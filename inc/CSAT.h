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
	CSAT(){}
	virtual ~CSAT(){}

	virtual bool GetKrw(int size,double *sw,double *krw) = 0 ;
	virtual bool GetKro(int size,double *so,double *kro) = 0 ;
	virtual double GetKrwPDSw(int size,double *sw,double *dkrw_dsw) = 0 ;
	virtual bool GetKroPDSw(int size,double *sw, double *dkrow_dsw) = 0 ;
	virtual double GetKroPDSg(double sg) = 0;
	virtual bool GetKroPDSg(int size, double *sg, double *dkrog_dsg) = 0 ;
	virtual bool GetKrgPDSg(int size,double *sg, double *dkrg_dsg) = 0 ;
	virtual bool GetKrgPDSw(int size,double *sw,double* dkrg_dsw) = 0 ;
	virtual bool GetKrg(int size, double *sg, double *krg) = 0 ;

	virtual void SetSWOF(int _nr_sat_tab_, char *satFile){}
	virtual	void SetSGOF(int _nr_sat_tab_,char *satFile){}
	virtual	bool CreatePDTable(int n,double *tab_y, double* tab_x, double *tab_dy_dx){
		return false;
	}

	void SetSoEndPoints(double sorg){sorg_ = sorg;}
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

protected:
	double swcr_, swco_, sw_max_;
	double sorg_;
	double sgcr_, sgco_;
};

#endif
