#include "CSTDProdWell.h"
#include "CState.h"
#include <iostream>
using namespace::std;

CSTDProdWell::CSTDProdWell(string well_name, int block_index)
:CStandardWell(well_name, block_index)
{
	well_type_ = STDPROD;
}

CSTDProdWell::~CSTDProdWell()
{

}
bool CSTDProdWell::CalParameterJacobian(CState *state){
	if (state->GetHasGas() == 1 && state->GetHasOil() == 1 && state->GetHasWater() == 0){
		dmo_dp_ = state->get_dmo_dpo(n_);
		dmo_dsg_ = state->get_dmo_dsg(n_);
		dmg_dp_ = state->get_dmg_dpg(n_);
		dmg_dsg_ = state->get_dmg_dsg(n_);
	}
	else if (state->GetHasGas() == 0 && state->GetHasOil() == 1 && state->GetHasWater() == 1){
		dmo_dp_ = state->get_dmo_dpo(n_);
		dmo_dsw_ = state->get_dmo_dsw(n_);
		dmw_dp_ = state->get_dmw_dpw(n_);
		dmw_dsw_ = state->get_dmw_dsw(n_);
		CalDqoDp(state);
		CalDqoDsw(state);
		CalDqwDp(state);
		CalDqwDsw(state);
		//
	}
	else{

	}
	return true;
}


bool CSTDProdWell::CalParameterResidual(CState *state){
//	cout << "Calculate Prod Well Parameters...." <<" gas:"<<state->GetHasGas()<<"water: "<<state->GetHasWater() << "oil: " <<state->GetHasOil();
	if (state->GetHasGas() == 1 && state->GetHasOil() == 1 && state->GetHasWater() == 0){
		po_ = state->get_po()[n_];
		pg_ = po_;
		mo_ = state->get_mo(n_);
		mg_ = state->get_mg(n_);

	}
	else if (state->GetHasGas() == 0 && state->GetHasOil() == 1 && state->GetHasWater() == 1){
		po_ = state->get_po()[n_];
		pw_ = po_;
		mo_ = state->get_mo(n_);
		mw_ = state->get_mw(n_);

		//
		CalORAT(state);
		CalBHP(state);
		CalWRAT(state);
		//

	}
	else{

	}
	return true;
}
bool CSTDProdWell::CalBHP(CState *state){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		BHP_ = TL_BHP_;
		break;
	case CStandardWell::CORAT:
		BHP_ = po_ - TL_ORAT_ / (mo_*WI_);
		break;
	case CStandardWell::CWRAT:
		BHP_ = pw_ - TL_WRAT_ / (mw_*WI_);
		break;
	case CStandardWell::CLRAT:
		BHP_ = (mo_*po_ + mw_*pw_ + mg_*pg_ + rso_*mo_*po_ - TL_LRAT_ / WI_) / (mo_ + mw_ + mg_ + rso_*mo_);
		break;
	case CStandardWell::CGRAT:
		BHP_ = (rso_*mo_*po_ + mg_*pg_ - TL_GRAT_ / WI_) / (rso_*mo_ + mg_);
		break;
	default:
		break;
	}
	return true;
}
bool CSTDProdWell::CalORAT(CState *state){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		ORAT_ = WI_*mo_*(po_ - TL_BHP_);
		break;
	case CStandardWell::CORAT:
		ORAT_ = TL_ORAT_;
		break;
	case CStandardWell::CWRAT:
		ORAT_ = (mo_ / mw_)*TL_WRAT_;
		break;
	case CStandardWell::CLRAT:
		ORAT_ = (mo_ / (mo_ + mw_ + mg_ + rso_*mo_))*TL_LRAT_;
		break;
	case CStandardWell::CGRAT:
		ORAT_ = (mo_ / (rso_*mo_ + mg_))*TL_GRAT_;
		break;
	default:
		break;
	}
	return true;
}
bool CSTDProdWell::CalWRAT(CState *state){
	int temp;
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		WRAT_ = WI_*mw_*(pw_ - TL_BHP_);
		break;
	case CStandardWell::CORAT:
//		cout << "mw:"<<mw_<<"mo:"<<mo_<<"TORAT"<<TL_ORAT_<<endl;
		WRAT_ = (mw_ / mo_)*TL_ORAT_;
//		cin >> temp ;
		break;
	case CStandardWell::CWRAT:
		WRAT_ = TL_WRAT_;
		break;
	case CStandardWell::CLRAT:
		WRAT_ = (mw_ / (mo_ + mw_ + mg_ + rso_*mo_))*TL_LRAT_;
		break;
	case CStandardWell::CGRAT:
		WRAT_ = (mw_ / (mo_ + rso_*mo_))*TL_GRAT_;
		break;
	default:
		break;
	}
	return true;
}
bool CSTDProdWell::CalGRAT(CState *state){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		GRAT_ = WI_*(mg_*(pg_ - TL_BHP_)+rso_*mo_*(po_-TL_BHP_));
		break;
	case CStandardWell::CORAT:
		GRAT_ = (rso_*mo_ + mg_)*TL_ORAT_ / mo_;
		break;
	case CStandardWell::CWRAT:
		GRAT_ = (rso_*mo_ + mg_)*TL_WRAT_ / mw_;
		break;
	case CStandardWell::CLRAT:
		GRAT_ = ((rso_*mo_ + mg_) / (mo_ + mw_ + mg_ + rso_*mo_))*TL_LRAT_;
		break;
	case CStandardWell::CGRAT:
		GRAT_ = TL_GRAT_;
		break;
	default:
		break;
	}
	return true;
}
bool CSTDProdWell::CalLRAT(CState *state){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		LRAT_ = WI_*(mg_*(pg_ - TL_BHP_) + (1 + rso_)*mo_*(po_ - TL_BHP_)+mw_*(pw_ - TL_BHP_));
		break;
	case CStandardWell::CORAT:
		LRAT_ = (rso_*mo_ + mg_ + mo_ + mw_)*TL_ORAT_ / mo_;
		break;
	case CStandardWell::CWRAT:
		LRAT_ = (rso_*mo_ + mg_ + mo_ + mw_)*TL_WRAT_ / mw_;
		break;
	case CStandardWell::CLRAT:
		LRAT_ = TL_LRAT_;
		break;
	case CStandardWell::CGRAT:
		LRAT_ = (rso_*mo_ + mg_ + mo_ + mw_)*TL_WRAT_ / (rso_*mo_ + mg_);
		break;
	default:
		break;
	}
	return true;
}

bool CSTDProdWell::CalDqoDp(CState *State){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		dqo_dp_ = WI_*(dmo_dp_*(po_ - BHP_) + mo_);
		break;
	case CStandardWell::CORAT:
		dqo_dp_ = 0;
		break;
	case CStandardWell::CWRAT:
		// To be completed
		break;
	case CStandardWell::CLRAT:
		// To be completed
		break;
	case CStandardWell::CGRAT:
		// To be completed
		break;
	default:
		break;
	}
	return true;
}
bool CSTDProdWell::CalDqwDp(CState *State){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		dqw_dp_ = WI_*(dmw_dp_*(pw_ - BHP_) + mw_);
		break;
	case CStandardWell::CORAT:
		dqw_dp_ = ORAT_*((mo_*dmw_dp_-mw_*dmo_dp_)/(mo_*mo_));
		break;
	case CStandardWell::CWRAT:
		dqw_dp_ = 0;
		break;
	case CStandardWell::CLRAT:
		// To be completed
		break;
	case CStandardWell::CGRAT:
		// To be completed
		break;
	default:
		break;
	}
	return true;
}
bool CSTDProdWell::CalDqgDp(CState *State){
	return true;
}
bool CSTDProdWell::CalDqoDsw(CState *State){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		dqo_dsw_ = WI_*dmo_dsw_*(po_ - BHP_);
		break;
	case CStandardWell::CORAT:
		dqo_dsw_ = 0;
		break;
	case CStandardWell::CWRAT:
		// To be completed
		break;
	case CStandardWell::CLRAT:
		// To be completed
		break;
	case CStandardWell::CGRAT:
		// To be completed
		break;
	default:
		break;
	}
	return true;
}
bool CSTDProdWell::CalDqwDsw(CState *State){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		dqw_dsw_ = WI_*dmw_dsw_*(pw_ - BHP_);
		break;
	case CStandardWell::CORAT:
		dqw_dsw_ = ORAT_*((mo_*dmw_dsw_ - mw_*dmo_dsw_) / (mo_*mo_));
		break;
	case CStandardWell::CWRAT:
		// To be completed
		break;
	case CStandardWell::CLRAT:
		// To be completed
		break;
	case CStandardWell::CGRAT:
		// To be completed
		break;
	default:
		break;
	}
	return true;
}

bool CSTDProdWell::CheckLimits(){
	switch (ctrl_mode_)
	{
	case CStandardWell::CBHP:
		if (ORAT_>TL_ORAT_){
			ORAT_ = TL_ORAT_;
			ctrl_mode_ = CORAT;
			return false;
		}
		break;
	case CStandardWell::CORAT:
		if (BHP_ < TL_BHP_){
			BHP_ = TL_BHP_;
			ctrl_mode_ = CBHP;
			return false;
		}
		break;
	case CStandardWell::CWRAT:
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
