#ifndef __CUnstructGrid_H__
#define __CUnstructGrid_H__

#include "CGrid.h"
#include <fstream>
#include <string>

using namespace::std;

class CUnstructGrid : public CGrid{
public:
		CUnstructGrid(int ncell) :CGrid(ncell){

		}
	  bool InputKx(double kx){return false;}
	  bool InputKx(double *kx){return false;}
	  bool InputKx(char *file_name){return false;}
	  bool InputKy(double ky){return false;}
	  bool InputKy(double *ky){return false;}
	  bool InputKy(char *file_name){return false;}
	  bool InputKz(double kz){return false;}
	  bool InputKz(double *kz){return false;}
	  bool InputKz(char *file_name){return false;}
	  int GetnX(){return 0;}
	  int GetnY(){return 0;}
	  int GetnZ(){return 0;}
	  bool InputDx(double dx){return false;}
	  bool InputDx(double *dx){return false;}
	  bool InputDy(double dy){return false;}
	  bool InputDy(double *dy){return false;}
	  bool InputDz(double dz){return false;}
	  bool InputDz(double *dz){return false;}
	  double GetDX(int n){return -1.0;}
	  double GetDY(int n){ return -1.0; }
	  double GetDZ(int n){ return -1.0; }
	  double GetKx(int n){ return -1.0; }
	  double GetKy(int n){ return -1.0; }
	  double GetKz(int n){ return -1.0; }

	  bool CalVolume(){ return false; }
	  bool CalDepth(){ return false; }
	  bool GenerateConnList(){ return false; }

	  int GetIndex(int i, int j, int k){ return 0; }
	  bool InputTops(double tops){ return false; }
};

#endif