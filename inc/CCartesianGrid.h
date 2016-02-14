#ifndef __CCartesianGrid_H__
#define __CCartesianGrid_H__

#include "CStructGrid.h"
#include <stdio.h>
//class CPVT;
//class CSAT;

class CCartGrid : public CStructGrid{
public:
	CCartGrid(int nx, int ny, int nz) :
		CStructGrid(nx, ny, nz){
		tops_ = new double[ncell_];
		memset(tops_, 0, ncell_*sizeof(double));
	}

	bool CalVolume(){
		for (int k = 0; k < ncell_; k++){
			v_[k] = dx_[k] * dy_[k] * dz_[k];
		}
		return true;
	}

	bool CalDepth(){
		int m = 0, n = 0;
		int topsIndex = 0;
		for (int k = 0; k < nz_; k++){
			for (int j = 0; j < ny_; j++){
				for (int i = 0; i < nx_; i++){
					if (k == 0){
						m = GetIndex(i, j, k);
						topsIndex = ny_*j + i;
						d_[m] = tops_[topsIndex] + 0.5*dz_[m];
					}
					else{
						n = GetIndex(i, j, k - 1);
						m = GetIndex(i, j, k);
						d_[m] = (d_[n] + 0.5*dz_[n] + 0.5*dz_[m]);
					}
					if (d_[m] - dz_[m] / 2 < dmin_){
						dmin_ = d_[m] - dz_[m] / 2;
					}
					if (d_[m] + dz_[m]>dmax_){
						dmax_ = d_[m] + dz_[m] / 2;
					}
				}
			}
		}
		return true;
	}

	bool GenerateConnList(){
		nconns_ = (nx_ - 1)*ny_*nz_ + nx_*(ny_ - 1)*nz_ + nx_*ny_*(nz_ - 1);
        conn_list_ = new ConnList(nconns_);
		int l = 0, r = 0;
		double trans = 0.0;
		int count = 0;
		// X direction connections
		for (int k = 0; k < nz_; k++){
			for (int j = 0; j < ny_; j++){
				for (int i = 0; i < nx_ - 1; i++){
					l = GetIndex(i, j, k);
					r = GetIndex(i + 1, j, k);
					trans = 2 * alpha *dy_[l] * dz_[l] / ((dx_[l] / kx_[l]) + (dx_[r] / kx_[r]));
                    conn_list_->InsertConnection(count,l,r,trans);
					count++;
				}
			}
		}
		// Y direction connections
		for (int k = 0; k < nz_; k++){
			for (int i = 0; i < nx_; i++){
				for (int j = 0; j < ny_ - 1; j++){
					l = GetIndex(i, j, k);
					r = GetIndex(i, j + 1, k);
					trans = 2 * alpha *dx_[l] * dz_[l] / ((dy_[l] / ky_[l]) + (dy_[r] / ky_[r]));
                    conn_list_->InsertConnection(count,l,r,trans);
					count++;
				}
			}
		}
		// Z direction connections
		for (int j = 0; j < ny_; j++){
			for (int i = 0; i < nx_; i++){
				for (int k = 0; k < nz_ - 1; k++){
					l = GetIndex(i, j, k);
					r = GetIndex(i, j, k + 1);
					trans = 2 * alpha *dx_[l] * dy_[l] / ((dz_[l] / kz_[l]) + (dz_[r] / kz_[r]));
                    conn_list_->InsertConnection(count,l,r,trans);
					count++;
				}
			}
		}

		return true;
	}

	int GetIndex(int i, int j, int k){
		int res = nx_*ny_*k + nx_*j + i;
		return res;
    }

	bool InputTops(double tops){
		for (int i = 0; i < ncell_; i++)
			tops_[i] = tops;
		return true;
    }
	~CCartGrid(){
		delete[] tops_;
	}
private:
	//Grid Variables
	double* tops_;
};





#endif
