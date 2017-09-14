#ifndef __CSimCtrl_H__
#define __CSimCtrl_H__


#include "CCartesianGrid.h"
#include "CPVT.h"
#include "CSAT.h"
#include "CSAT_TABLE.h"
#include "CSAT_COREY.hpp"
#include "CSolver.h"
#include "CSchedule.h"
#include "CSTDProdWell.h"
#include "CSTDWInjWell.h"
#include "CHistoryMatching.hpp"
#include <vector>

using namespace std;

class CState;

class SimCtrl{
public:
	SimCtrl();
	~SimCtrl();
public:
	
	enum GRIDTYPE{CARTESIAN,CORNERPOINT,UNSTRUCTURED};
	enum MODELTYPE{SBLACKOIL,COMPONENT};
	enum PHASETYPE{OWG,OW,OG};
	enum IMPLICTLEVEL{FIM,IMPES,AIM};
	enum LINEARSOLVER{PARDISO,GMRES};
	enum PRECONDITIONER{NA,TRUECPR,ILU0};

public:
	
	bool set_grid_type(GRIDTYPE grid_type);
	bool set_model_type(MODELTYPE model_type);
	bool set_phase_type(PHASETYPE phase_type);
	bool set_implict_level(IMPLICTLEVEL implict_level);
	bool set_linear_solver(LINEARSOLVER linear_solver);
	bool set_preconditioner(PRECONDITIONER preconditioner);

public:
	bool InitializeRun(char* input_file);
	bool InitializeRun();
	bool InitializeGrid();
	bool InitializeGrid(int ncell);
	bool CheckInputValidity();
	bool InitializePVT();
	bool InitializeSAT();
	bool InitializeSCH();
	bool InitializeState();
	bool InitializeSolver();
	bool InitializeWell();

public:
	bool SetInitPres(double datumDepth, double pDatum);
	bool SetInitPres(double *po);
	bool SetInitSat(double Sw, double Sg);
	bool SetInitSat(double *Sw, double *Sg);
public:
	bool RunSim();
public:
	bool OutputResult();
public:
	bool Display();
public:
	// Options

	GRIDTYPE grid_type_;
	MODELTYPE model_type_;
	PHASETYPE phase_type_;
	IMPLICTLEVEL implict_level_;
	LINEARSOLVER linear_solver_;
	PRECONDITIONER preconditioner_;
	// Newton-Raphson Options
	int max_newton_iters_;
	double min_err;

public:
	// Data
	CGrid* grid_;
//	CCartGrid* grid_;
	CPVT* pvt_;
	CSAT* sat_;
	CState* state_;
	CSchedule* sch_;
	CSolver* solver_;
	vector<CStandardWell*> std_well_;
	CHistoryMatching* hm_;
	//	CProdWell* prod_;

	//Newton Iteration performance
	vector<int> tstep_;
	vector<double> time_;
	vector<int> num_newton_iter_;
	vector<double> cfl_;
	vector<double> dt_;

public:
	int display_level_;
};

#endif
