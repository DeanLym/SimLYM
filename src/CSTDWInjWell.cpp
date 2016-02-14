#include "CSTDWInjWell.h"
#include "CState.h"
#include <iostream>
using namespace::std;


CSTDWInjWell::CSTDWInjWell(string well_name, int block_index)
:CStandardWell(well_name, block_index)
{
	well_type_ = STDWINJ;
}

CSTDWInjWell::~CSTDWInjWell()
{

}
bool CSTDWInjWell::CalParameterJacobian(CState *state){
	if (state->GetHasGas() == 1 && state->GetHasOil() == 1 && state->GetHasWater() == 0){

	}
	else if (state->GetHasGas() == 0 && state->GetHasOil() == 1 && state->GetHasWater() == 1){
		dmo_dp_ = state->get_dmo_dpo(n_);
		dmo_dsw_ = state->get_dmo_dsw(n_);
		dmw_dp_ = state->get_dmw_dpw(n_);
		dmw_dsw_ = state->get_dmw_dsw(n_);
		//
		CalDqwDp(state);
		CalDqwDsw(state);
		//
	}
	else{

	}
	return true;
}
bool CSTDWInjWell::CalParameterResidual(CState *state){
//	cout << "Calculating Injector Parameters" << endl;
	if (state->GetHasGas() == 1 && state->GetHasOil() == 1 && state->GetHasWater() == 0){

	}
	else if (state->GetHasGas() == 0 && state->GetHasOil() == 1 && state->GetHasWater() == 1){
		pw_ = state->get_pw()[n_];
		po_ = pw_;
		mo_ = state->get_mo(n_);
		mw_ = state->get_mw(n_);

		//
		CalBHP(state);
		CalWRAT(state);

	}
	else{

	}
	return true;
}
bool CSTDWInjWell::CalBHP(CState *state){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		BHP_ = TL_BHP_;
		break;
	case CStandardWell::CORAT:
		return false;
		break;
	case CStandardWell::CWRAT:
		BHP_ = pw_ - TL_WRAT_ / (WI_*(mo_ + mw_ + mg_));
		break;
	case CStandardWell::CLRAT:
		return false;
		break;
	case CStandardWell::CGRAT:
		return false;
		break;
	default:
		break;
	}
	return true;

}
bool CSTDWInjWell::CalORAT(CState *state){
	return true;
}
bool CSTDWInjWell::CalWRAT(CState *state){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		WRAT_ = WI_*(mo_ + mw_ + mg_)*(pw_ - TL_BHP_);
//		WRAT_ = WI_*mw_*(pw_ - TL_BHP_);

		break;
	case CStandardWell::CORAT:
		return false;
		break;
	case CStandardWell::CWRAT:
		WRAT_ = TL_WRAT_;
		break;
	case CStandardWell::CLRAT:
		return false;
		break;
	case CStandardWell::CGRAT:
		return false;
		break;
	default:
		break;
	}
	return true;
}
bool CSTDWInjWell::CalGRAT(CState *state){
	return false;
}
bool CSTDWInjWell::CalLRAT(CState *state){
	return false;
}
bool CSTDWInjWell::CalDqoDp(CState *State){
	return false;
}
bool CSTDWInjWell::CalDqwDp(CState *State){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		dqw_dp_ = WI_*((dmo_dp_ + dmw_dp_)*(pw_ - TL_BHP_) + (mo_ + mw_));
		break;
	case CStandardWell::CORAT:
		return false;
		break;
	case CStandardWell::CWRAT:
		dqw_dp_ = 0;
		break;
	case CStandardWell::CLRAT:
		return false;
		break;
	case CStandardWell::CGRAT:
		return false;
		break;
	default:
		break;
	}
	return true;
}
bool CSTDWInjWell::CalDqgDp(CState *State){
	return false;
}
bool CSTDWInjWell::CalDqoDsw(CState *State){
	return false;
}
bool CSTDWInjWell::CalDqwDsw(CState *State){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		dqw_dsw_ = WI_*(dmo_dsw_ + dmw_dsw_)*(pw_ - TL_BHP_);
		break;
	case CStandardWell::CORAT:
		return false;
		break;
	case CStandardWell::CWRAT:
		dqw_dsw_ = 0;
		break;
	case CStandardWell::CLRAT:
		return false;
		break;
	case CStandardWell::CGRAT:
		return false;
		break;
	default:
		break;
	}
	return true;
}

bool CSTDWInjWell::CheckLimits(){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		if (WRAT_>TL_WRAT_){
			WRAT_ = TL_WRAT_;
			ctrl_mode_ = CWRAT;
			return false;
		}
		break;
	case CStandardWell::CORAT:
		break;
	case CStandardWell::CWRAT:
		if (BHP_ > TL_BHP_){
			BHP_ = TL_BHP_;
			ctrl_mode_ = CBHP;
			return false;
		}
		break;
	case CStandardWell::CLRAT:
		break;
	case CStandardWell::CGRAT:
		break;
	default:
		break;
	}
	return true;
}
