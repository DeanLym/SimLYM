#ifndef __CPVT_H__
#define __CPVT_H__

#include <stdio.h>
#include <string>
#include <math.h>
#include <fstream>

using namespace std;
using std::string;


class CPVT{
public:

	CPVT()
:po_tab_(NULL),
 bo_tab_(NULL),
 rso_tab_(NULL),
 vo_tab_(NULL),
 pg_tab_(NULL),
 vg_tab_(NULL),
 bg_tab_(NULL),
 dbo_dpo_tab_(NULL),
 dvo_dpo_tab_(NULL),
 drso_dpo_tab_(NULL),
 dbg_dpg_tab_(NULL),
 dvg_dpg_tab_(NULL),
 cv_(0.0)
{

}

	~CPVT(){
		delete []po_tab_; po_tab_ = NULL;
		delete []bo_tab_; bo_tab_ = NULL;
		delete []vo_tab_; vo_tab_ = NULL;
		delete []rso_tab_; rso_tab_ = NULL;
		delete []pg_tab_; pg_tab_ = NULL;
		delete []vg_tab_; vg_tab_ = NULL;
		delete []bg_tab_; bg_tab_ = NULL;
		delete []dbo_dpo_tab_; dbo_dpo_tab_ = NULL;
		delete []dvo_dpo_tab_; dvo_dpo_tab_ = NULL;
		delete []drso_dpo_tab_; drso_dpo_tab_ = NULL;
		delete []dbg_dpg_tab_; dbg_dpg_tab_ = NULL;
		delete []dvg_dpg_tab_; dvg_dpg_tab_ = NULL;

	}
	void set_buble_point_pressure(double buble_point_pressure){
		buble_point_pressure_ = buble_point_pressure;
	}
	int SetDensity(double dO,double dW,double dG){
		denO_std = dO;
		denW_std = dW;
		denG_std = dG;
		return 1;
	}
	bool GetBo(int size,double *po,double* bo){
		//		double co = 1.2e-5;
		//		for(int i = 0; i < size; i++){
		//			bo[i] = exp(-co*(po[i] - p_ref_));
		//        }
		// Table lookup
		for (int i=0;i < size;i++){
			bo[i] = PiecewiseLinearInter(nr_pvt_oil_,po_tab_,bo_tab_,po[i]);
		}
		return true;
	}

	double GetBo(double po){
		double bo = PiecewiseLinearInter(nr_pvt_oil_,po_tab_,bo_tab_,po);
		return bo;
	}

	bool GetBw(int size, double* pw,double* bw){
		//		double cw = 5e-7;
		//		for(int i = 0; i < size; i++){
		//			bw[i] = exp(-cw*(pw[i] - p_ref_));
		//		}
		double x;
		for (int i=0;i < size;i++){
			x = cw_*(pw[i]-pw_ref_);
			bw[i] = bw_ref_/(1+x+x*x/2);
		}
		return true;
	}

	bool GetVisO(int size, double* po,double* vo){
		for (int i=0;i < size;i++){
			vo[i] = PiecewiseLinearInter(nr_pvt_oil_,po_tab_,vo_tab_,po[i]);
		}
		return true;
	}

	bool GetVisW(int size, double* pw,double* vw){
		double y;
		if(cv_ == 0.0)
			for(int i=0;i<size;i++)
				vw[i] = vw_ref_;
		else
			for(int i=0;i<size;i++){
				y = -cv_*(pw[i]-pw_ref_);
				vw[i] = vw_ref_/(1+y+y*y/2);
			}
		return true;
	}

	bool GetBoPD(int size, double* po,double* dbo_dpo_){
		for (int i=0;i < size;i++){
			dbo_dpo_[i] = PiecewiseLinearInter(nr_pvt_oil_,po_tab_,dbo_dpo_tab_,po[i]);
		}
		return true;
	}

	bool GetVisOPD(int size, double* po,double* dvo_dpo_){
		for (int i=0;i < size;i++){
			dvo_dpo_[i] = PiecewiseLinearInter(nr_pvt_oil_,po_tab_,dvo_dpo_tab_,po[i]);
		}
		return true;
	}

	bool GetBwPD(int size, double* pw,double* dbw_dpw_){
		double x, deno;
		for (int i=0;i < size;i++){
			x = cw_*(pw[i]-pw_ref_);
			deno = 1+x+x*x/2;
			dbw_dpw_[i] = -bw_ref_*cw_*(1+x)/(deno*deno);
		}
		return true;
	}

	bool GetVisWPD(int size, double* pw,double* dvw_dpw_){
		double y, deno;
		if(cv_ == 0.0)
			for(int i=0;i<size;i++)
				dvw_dpw_[i] = 0.0;
		else
			for (int i=0;i < size;i++){
				y = -cv_*(pw[i]-pw_ref_);
				deno = 1+y+y*y/2;
				dvw_dpw_[i] = vw_ref_*cv_*(1+y)/(deno*deno);
			}
		return true;
	}

	bool GetDenO(int size, double* bo, double* rso, double* deno){
		for (int i = 0; i<size; i++){
			deno[i] = (denO_std + rso[i] * denG_std) / bo[i];
		}
		return true;
	}
	double GetDenO(double po){
		double deno;
		double bo;
		bo = GetBo(po);
		deno = denO_std / bo;
		return deno;
	}
	bool GetDenO(int size, double* bo, double* deno){
		for (int i = 0; i<size; i++){
			deno[i] = denO_std / bo[i];
		}
		return true;
	}
	bool GetDenW(int size, double* bw, double* denw){
		for (int i = 0; i<size; i++){
			denw[i] = denW_std / bw[i];
		}
		return true;
	}
	bool GetDenG(int size, double *bg, double* deng){
		for (int i = 0; i<size; i++){
			deng[i] = denG_std / bg[i];
		}
		return true;
	}

	double GetDenO_Std(){ return denO_std; }
	double GetDenW_Std(){ return denW_std; }
	double GetDenG_Std(){ return denG_std; }

	bool GetBg(int size, double* pg, double* bg){
		for (int i = 0; i < size; i++){
			bg[i] = 3592.7*pow(pg[i], -1.0226) / 178.1;
		}
		return true;
	}
	bool GetBgPD(int size, double* pg, double* dbg_dpg_){
		for (int i = 0; i<size; i++){
			dbg_dpg_[i] = (-1.0226*3592.7 / 178.1)*pow(pg[i], -2.0226);
		}
		return true;
	}
	bool GetVisG(int size, double* pg, double* vg){
		double temp = 0.0;
		for (int i = 0; i<size; i++){
			temp = pg[i];
			vg[i] = (3e-10)*temp*temp + (1e-6)*temp + 0.0133;
		}
		return true;
	}
	bool GetVisGPD(int size, double* pg, double* dvg_dpg_){
		for (int i = 0; i<size; i++){
			dvg_dpg_[i] = (6e-10)*pg[i] + (1e-6);
		}
		return true;
	}
	bool GetRso(int size, double* po, double* rso){
		for (int i = 0; i < size; i++){
			if (po[i] > buble_point_pressure_){
				rso[i] = 178.1;
			}
			else{
				rso[i] = 178.1*0.0000476*pow(po[i], 1.206511);
			}
		}
		//		res = 178.1*PiecewiseLinearInter(nr_pvt_oil_, po_tab_, rso_tab_, po);
		return true;
	}

	bool GetRsoPD(int size, double* po, double* drso_dpo){
		double pBub = 3824.321712;
		for (int i = 0; i<size; i++){
			if (po[i] > pBub){
				drso_dpo[i] = 0.0;
			}
			else{
				drso_dpo[i] = 178.1*0.0000476*1.206511*pow(po[i], 0.206511);
			}
		}
		return true;
	}

	void set_pvdo_table(int num_node_in_pvdo_table, char *pvdo_file){
		nr_pvt_oil_ = num_node_in_pvdo_table;
		po_tab_ = new double[nr_pvt_oil_];
		vo_tab_ = new double[nr_pvt_oil_];
		bo_tab_ = new double[nr_pvt_oil_];
		dbo_dpo_tab_ = new double[nr_pvt_oil_];
		dvo_dpo_tab_ = new double[nr_pvt_oil_];

		fstream fp;
		fp.open(pvdo_file, ios::in);

		string temp;
		fp >> temp;

		for (int i = 0; i < nr_pvt_oil_; i++){
			fp >> po_tab_[i];
			fp >> bo_tab_[i];
			fp >> vo_tab_[i];
		}
		fp.close();

		CreatePDTable(nr_pvt_oil_,bo_tab_,po_tab_,dbo_dpo_tab_);
		CreatePDTable(nr_pvt_oil_,vo_tab_,po_tab_,dvo_dpo_tab_);
	}

	void set_pvto_table(int num_node_in_pvto_table, char *pvto_file){
		nr_pvt_oil_ = num_node_in_pvto_table;
		po_tab_ = new double[nr_pvt_oil_];
		vo_tab_ = new double[nr_pvt_oil_];
		bo_tab_ = new double[nr_pvt_oil_];
		rso_tab_ = new double[nr_pvt_oil_];
		dbo_dpo_tab_ = new double[nr_pvt_oil_];
		dvo_dpo_tab_ = new double[nr_pvt_oil_];
		drso_dpo_tab_ = new double[nr_pvt_oil_];

		fstream fp;
		fp.open(pvto_file, ios::in);

		string temp;
		fp >> temp;

		for (int i = 0; i < nr_pvt_oil_; i++){
			fp >> rso_tab_[i];
			fp >> po_tab_[i];
			fp >> bo_tab_[i];
			fp >> vo_tab_[i];
		}
		fp.close();

		CreatePDTable(nr_pvt_oil_,bo_tab_,po_tab_,dbo_dpo_tab_);
		CreatePDTable(nr_pvt_oil_,vo_tab_,po_tab_,dvo_dpo_tab_);
		CreatePDTable(nr_pvt_oil_,rso_tab_,po_tab_,drso_dpo_tab_);
	}

	void set_pvtw_table(char* pvtw_file){
		fstream fp;
		fp.open(pvtw_file, ios::in);
		string temp;
		fp >> temp;

//		for(int i=0;i<5;i++)
			fp >> pw_ref_ >> bw_ref_ >> cw_ >> vw_ref_ >> cv_;

		fp.close();

	}
	void set_pvdg_table(int num_node_in_pvdg_table, char *pvdg_file){
		nr_pvt_gas_ = num_node_in_pvdg_table;
		pg_tab_ = new double[nr_pvt_gas_];
		vg_tab_ = new double[nr_pvt_gas_];
		bg_tab_ = new double[nr_pvt_gas_];
		dbg_dpg_tab_ = new double[nr_pvt_gas_];
		dvg_dpg_tab_ = new double[nr_pvt_gas_];

		fstream fp;
		fp.open(pvdg_file, ios::in);

		string temp;
		fp >> temp;

		for (int i = 0; i < nr_pvt_gas_; i++){
			fp >> pg_tab_[i];
			fp >> bg_tab_[i];
			fp >> vg_tab_[i];
		}
		fp.close();

		CreatePDTable(nr_pvt_gas_,bg_tab_,pg_tab_,dbg_dpg_tab_);
		CreatePDTable(nr_pvt_gas_,vg_tab_,pg_tab_,dvg_dpg_tab_);
	}
	bool CreatePDTable(int n,double *tab_y, double* tab_x, double *tab_dy_dx){
		tab_dy_dx[0] = (tab_y[1] - tab_y[0])/(tab_x[1]-tab_x[0]);
		for(int i=1 ; i<n-1 ; i++){
			tab_dy_dx[i] = (tab_y[i+1] - tab_y[i-1])/(tab_x[i+1]-tab_x[i-1]);
		}
		tab_dy_dx[n-1] =
				(tab_y[n -1 ] - tab_y[n -2])
				/(tab_x[n -1 ]-tab_x[n -2]);
		return true;
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
private:

	//	int numOPressNodes,numWPressNodes,numGasPressNodes;
	double denO_std, denW_std, denG_std;
	//	double *pPo, *pPw;
	//	double *pBo, *pBw, *pBg;
	//	double *pRso;
	//	double *pVisO, *pVisW, *pVisG;
	int nr_pvt_oil_;
	double *po_tab_, *vo_tab_;
	double *bo_tab_, *rso_tab_;
	double *dbo_dpo_tab_, *dvo_dpo_tab_, *drso_dpo_tab_;
	double *dbg_dpg_tab_, *dvg_dpg_tab_;
	double buble_point_pressure_;
	double pw_ref_,bw_ref_,cw_, vw_ref_, cv_;

	const double p_ref_ = 14.7;

	int nr_pvt_gas_;

	double *pg_tab_, *vg_tab_, *bg_tab_;

};

#endif
