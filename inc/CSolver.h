#ifndef __CSOLVER_H__
#define __CSOLVER_H__
#include"CMatrixBlock.h"
class MatrixBlock;
class CState;

class CSolver{
public:
	CSolver(CState* state){
		state_ = state;
	}
	virtual const double* Solve() const= 0;
	virtual bool InputResidual(const double* resid) = 0;
	virtual bool InputJacobian() = 0;
	~CSolver(){

	}
private:
	CState* state_;
	
};


#endif