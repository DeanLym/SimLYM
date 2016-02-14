#include "CPardisoSolver.h"
#include "Interface_Pardiso.h"
#include "CState.h"
#include "CMatrixBlock.h"
#include <algorithm>

using namespace::std;


struct JACO2CRS{
	int n;
	int row;
	int col;
};

bool comp(const JACO2CRS& lhs, const JACO2CRS& rhs){
	if (lhs.row < rhs.row){
		return true;
	}
	else if (lhs.row == rhs.row){
		if (lhs.col < rhs.col){
			return true;
		}
		else{
			return false;
		}
	}
	else{
		return false;
	}
}


CPardisoSolver::CPardisoSolver(CState* state)
:CSolver(state)
{
	int num_phase = state->get_num_phase();
	int nconns = state->get_nconns();
	int ncells = state->get_n_act_cell();
	block_size_ = num_phase*num_phase;
	nnz_ = block_size_*(ncells + 2 * nconns);
	nrow_ = num_phase*ncells;
	int num_d = block_size_*ncells;
	int num_od = block_size_*nconns;
	int* n_ = new int[nnz_];
	pIA_ = new int[nrow_ + 1];
	pJA_ = new int[nnz_];
	pA_ = new double[nnz_];
	pB_ = new double[nrow_];
	index_ = new const double*[nnz_];
	const int* l = state->get_l();
	const int* r = state->get_r();
	const MatrixBlock* jaco_d = state->get_jaco_d();
	const MatrixBlock* jaco_ld = state->get_jaco_ld();
	const MatrixBlock* jaco_ud = state->get_jaco_ud();

	JACO2CRS* a = new JACO2CRS[nnz_];
//	JACO2CRS a[nnz];
	int count = 0;

	for (int i = 0; i < ncells; i++){
		for (int j = 0; j < num_phase; j++){
			for (int k = 0; k < num_phase; k++){
				a[count].n = count;
				a[count].row = num_phase*i + j;
				a[count].col = num_phase*i + k;
				count++;
			}
		}
	}

	for (int i = 0; i < nconns; i++){
		for (int j = 0; j < num_phase; j++){
			for (int k = 0; k < num_phase; k++){
				a[count].n = count;
				a[count].row = num_phase*l[i] + j;
				a[count].col = num_phase*r[i] + k;
				count++;
			}
		}
	}

	for (int i = 0; i < nconns; i++){
		for (int j = 0; j < num_phase; j++){
			for (int k = 0; k < num_phase; k++){
				a[count].n = count;
				a[count].row = num_phase*r[i] + j;
				a[count].col = num_phase*l[i] + k;
				count++;
			}
		}
	}
//	sort(begin(a), end(a), comp);
	sort(a, a+nnz_, comp);
	count = 0;
	for (int i = 0; i < nnz_; i++){
		if (i == 0){
			pIA_[count] = i + 1;
			count++;
		}
		else if (a[i].row>a[i - 1].row){
			pIA_[count] = i + 1;
			count++;
		}
		n_[i] = a[i].n;
		pJA_[i] = a[i].col + 1;
	}
	pIA_[count] = nnz_ + 1;


	int n, k, m;
	for (int i = 0; i < nnz_; i++){
		n = n_[i];
		if (n < num_d){  //Located in Jaco_d
			k = n / block_size_;
			m = n % block_size_;
			index_[i] = jaco_d[k].get_pointer(m);
		}
		else if (n < num_d + num_od){ // Located in Jaco_ud
			k = (n - num_d) / block_size_;
			m = (n - num_d) % block_size_;
			index_[i] = jaco_ud[k].get_pointer(m);
		}
		else{ //Located in Jaco_ld
			k = (n - num_d - num_od) / block_size_;
			m = (n - num_d - num_od) % block_size_;
			index_[i] = jaco_ld[k].get_pointer(m);
		}
	}

	int stop = 0;
	delete[] n_;
}

const double* CPardisoSolver::Solve() const{
	UGSOL::MKLPardiso PardisoSolver;
	const int idx = 1;
	const int mtype = 11;
	const int mnum = 1;
	const int cgs = 2;
	const int msglvl = 0;
//	cout << "Initialize Pardiso" << endl;
	PardisoSolver.Init_Pardiso(nrow_, pA_, pIA_, pJA_, idx, mtype, mnum, cgs, msglvl);
//	cout << endl<<"Pardiso Solve" << endl;
	PardisoSolver.Pardiso_Solve(nrow_, pA_, pIA_, pJA_, pB_);

	return pB_;
}
bool CPardisoSolver::InputResidual(const double* resid){
	for (int i = 0; i < nrow_; i++){
		pB_[i] = -resid[i];
	}
	return true;
}
bool  CPardisoSolver::InputJacobian(){
	const double* index;
	for (int i = 0; i < nnz_; i++){
		index = index_[i];
		pA_[i] = *index;
	}
	return true;
}

bool CPardisoSolver::Initialize(CState *state){




	return true;

}


CPardisoSolver::~CPardisoSolver(){
	delete[] pIA_;
	delete[] pJA_;
	delete[] pA_;
	delete[] pB_;
	delete[] index_;
}
