#ifndef __CSTDINJWELLSP_H__
#define __CSTDINJWELLSP_H__
#include "CStandardWell.h"
class CSTDWInjWell: public CStandardWell
{
public:
	enum INJPHASE{INJWATER,INJGAS};
	CSTDWInjWell(string well_name, int block_index);
	~CSTDWInjWell();
public:
	bool CalParameterJacobian(CState *state);
	bool CalParameterResidual(CState *state);
	bool CalBHP(CState *state);
	bool CalORAT(CState *state);
	bool CalWRAT(CState *state);
	bool CalGRAT(CState *state);
	bool CalLRAT(CState *state);
	bool CalDqoDp(CState *State);
	bool CalDqwDp(CState *State);
	bool CalDqgDp(CState *State);
	bool CalDqoDsw(CState *State);
	bool CalDqwDsw(CState *State);
public:
	bool CheckLimits();

	

private:




};




#endif