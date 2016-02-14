#ifndef __CStructGrid_H__
#define __CStructGrid_H__

#include "CGrid.h"
#include <fstream>
#include <string>

using namespace::std;

class CStructGrid : public CGrid{
public:
	CStructGrid(int nx, int ny, int nz) :
		CGrid(nx*ny*nz){
		nx_ = nx;
		ny_ = ny;
		nz_ = nz;

		kx_ = new double[ncell_];
		memset(kx_, 0, ncell_*sizeof(double));
		ky_ = new double[ncell_];
		memset(ky_, 0, ncell_*sizeof(double));
		kz_ = new double[ncell_];
		memset(kz_, 0, ncell_*sizeof(double));
		dx_ = new double[ncell_];
		memset(dx_, 0, ncell_*sizeof(double));
		dy_ = new double[ncell_];
		memset(dy_, 0, ncell_*sizeof(double));
		dz_ = new double[ncell_];
		memset(dz_, 0, ncell_*sizeof(double));
	}

	bool InputKx(double kx){
		for (int i = 0; i < ncell_; i++)
			kx_[i] = kx;
		return true;
	};
	bool InputKx(double *kx){
		for (int i = 0; i < ncell_; i++)
			kx_[i] = kx[i];
		return true;
	};

	bool InputKx(char *file_name){
		ifstream in;
		in.open(file_name);
		string temp;
		in >> temp;
		for (int i = 0; i < ncell_; i++){
			in >> kx_[i];
		}
		in.close();
		return true;
	}
	bool InputKy(double ky){
		for (int i = 0; i < ncell_; i++)
			ky_[i] = ky;
		return true;
	};
	bool InputKy(char *file_name){
		ifstream in;
		in.open(file_name);
		string temp;
		in >> temp;
		for (int i = 0; i < ncell_; i++){
			in >> ky_[i];
		}
		in.close();
		return true;
	}
	bool InputKy(double *ky){
		for (int i = 0; i < ncell_; i++)
			ky_[i] = ky[i];
		return true;
	};
	bool InputKz(double kz){
		for (int i = 0; i < ncell_; i++)
			kz_[i] = kz;
		return true;
	};
	bool InputKz(char *file_name){
		ifstream in;
		in.open(file_name);
		string temp;
		in >> temp;
		for (int i = 0; i < ncell_; i++){
			in >> kz_[i];
		}
		in.close();
		return true;
	}
	bool InputKz(double *kz){
		for (int i = 0; i < ncell_; i++)
			kz_[i] = kz[i];
		return true;
	};
	bool InputDx(double dx){
		for (int i = 0; i < ncell_; i++)
			dx_[i] = dx;
		return true;
	}
	bool InputDx(double *dx){
		for (int i = 0; i < ncell_; i++)
			dx_[i] = dx[i];
		return true;
	};

	bool InputDy(double dy){
		for (int i = 0; i < ncell_; i++)
			dy_[i] = dy;
		return true;
	}
	bool InputDy(double *dy){
		for (int i = 0; i < ncell_; i++)
			dy_[i] = dy[i];
		return true;
	};
	bool InputDz(double dz){
		for (int i = 0; i < ncell_; i++)
			dz_[i] = dz;
		return true;
	}
	bool InputDz(double *dz){
		for (int i = 0; i < ncell_; i++)
			dz_[i] = dz[i];
		return true;
	};
	int GetnX(){
		return nx_;
	};
	int GetnY(){
		return ny_;
	};
	int GetnZ(){
		return nz_;
	};
	double GetDX(int n){ return dx_[n]; }
	double GetDY(int n){ return dy_[n]; }
	double GetDZ(int n){ return dz_[n]; }
	double GetKx(int n){ return kx_[n]; }
	double GetKy(int n){ return ky_[n]; }
	double GetKz(int n){ return kz_[n]; }
	~CStructGrid(){
		delete[] kx_;
		delete[] ky_;
		delete[] kz_;
		delete[] dx_;
		delete[] dy_;
		delete[] dz_;
	}
protected:
	int nx_;
	int ny_;
	int nz_;
	double *kx_, *ky_, *kz_;
	double *dx_, *dy_, *dz_;
};












#endif