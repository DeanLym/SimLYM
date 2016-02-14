#ifndef __CSTDPRODWELL_H__
#define __CSTDPRODWELL_H__

#include "CStandardWell.h"

class CSTDProdWell :public CStandardWell
{
public:
	CSTDProdWell(string well_name, int block_index);
	~CSTDProdWell();
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
};




#endif