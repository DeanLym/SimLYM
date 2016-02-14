#include "CSolver.h"

class CPardisoSolver : public CSolver
{
public:
	CPardisoSolver(CState* state);
	const double* Solve() const;
	bool InputResidual(const double* resid);
	bool  InputJacobian();
	~CPardisoSolver();
	bool Initialize(CState *state);
private:
//	int* n_;
	const double **index_;
	int nnz_, nrow_, block_size_;

	int *pIA_, *pJA_;
	double *pA_;
	double *pB_;
//	double pA_[nnz];
//	double pB_[nrow];
//	int pIA_[nrow + 1];
//	int pJA_[nnz];
};
