/*
 * CSAT_COREY.hpp
 *
 *  Created on: Mar 6, 2016
 *      Author: yiminliu
 */

#ifndef CSAT_COREY_HPP_
#define CSAT_COREY_HPP_

#include "CSAT.h"
#include <stdio.h>
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <exception>
#include <stdexcept>

using namespace std;

class CSAT_COREY:public CSAT{
public:
	CSAT_COREY():
	swi_(0.0),
	sor_(0.0),
	aw_(1.0),
	ao_(1.0),
	krw0_(1.0),
	kro0_(1.0),
	ds_(1.0)
	{

	}
	CSAT_COREY(double swi,double sor,double aw,double ao,double krw0,double kro0){
		swi_ = swi;
		sor_ = sor;
		aw_ = aw;
		ao_ = ao;
		krw0_ = krw0;
		kro0_ = kro0;
		ds_ = 1 - swi_ - sor_;
	}
	~CSAT_COREY(){

	}

	bool GetKrw(int size,double *sw,double *krw){
		double sw_norm;
		for(int i=0; i<size; i++){
			if(sw[i] <= swi_)
				krw[i] = 0.0;
			else if(sw[i] >= 1-sor_)
				krw[i] = krw0_;
			else{
				sw_norm = (sw[i] - swi_)/(1-swi_-sor_);
				krw[i] = krw0_ * pow(sw_norm,aw_);
			}
		}
		return true;
	}
	bool GetKro(int size,double *so,double *kro){
		double so_norm;
		for(int i=0; i<size; i++){
			if(so[i] <= sor_)
				kro[i] = 0.0;
			else if(so[i] >= 1-swi_)
				kro[i] = kro0_;
			else{
				so_norm = (so[i] - sor_)/ds_;
				kro[i] = kro0_ * pow(so_norm,ao_);
			}
		}
		return true;
	}

	double GetKrwPDSw(int size,double *sw,double *dkrw_dsw){
		double sw_norm;
		for(int i=0; i<size; i++){
			if(sw[i] <= swi_ || sw[i]>= 1-sor_)
				dkrw_dsw[i] = 0.0;
			else{
				sw_norm = (sw[i] - swi_)/ds_;
				dkrw_dsw[i] = krw0_ *aw_ / ds_* pow(sw_norm,aw_-1);
			}
		}
		return true;
	}

	bool GetKroPDSw(int size,double *sw, double *dkrow_dsw){
		double so_norm;
		for(int i=0; i<size; i++){
			if(sw[i] <= swi_ || sw[i]>= 1-sor_)
				dkrow_dsw[i] = 0.0;
			else{
				so_norm = (1 - sw[i] - sor_)/ds_;
				dkrow_dsw[i] = -kro0_ *ao_ / ds_* pow(so_norm,ao_ - 1);
			}
		}
		return true;
	}

	double GetKroPDSg(double sg){
		throw runtime_error("Incomplete implementation called.");
		return 1.0;
	}

	bool GetKroPDSg(int size, double *sg, double *dkrog_dsg){
		throw runtime_error("Incomplete implementation called.");
		return true;
	}

	bool GetKrgPDSg(int size,double *sg, double *dkrg_dsg){
		throw runtime_error("Incomplete implementation called.");
		return true;
	}

	bool GetKrgPDSw(int size,double *sw,double* dkrg_dsw){
		throw runtime_error("Incomplete implementation called.");
		return true;
	}
	bool GetKrg(int size, double *sg, double *krg){
		throw runtime_error("Incomplete implementation called.");
		return true;
	}

public:
	double ds_;
	double swi_, sor_;
	double aw_, ao_;
	double krw0_, kro0_;

};




#endif /* CSAT_COREY_HPP_ */
