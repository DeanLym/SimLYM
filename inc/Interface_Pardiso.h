#ifndef Interface_Pardiso_h__
#define Interface_Pardiso_h__

#include  "mkl.h"

namespace UGSOL {

#define PSIZE 64

class MKLPardiso {
public:
	MKLPardiso();
    ~MKLPardiso();
	void Clear_Pardiso();
	int Pardiso_Solve (int n, double* a, int* ia, int* ja, double* b);
	int Init_Pardiso  (int n, double* a, int* ia, int* ja, const int idx_base, 
                       const int _mtype, const int _mnum, const int _cgs, const int _msglvl);
	int Pardiso_Factor(int n, double* a, int* ia, int* ja);	
	int BackSubstitute(int n, double* a, int* ia, int* ja, double* b);

	double* solution;

private:
	int msglvl;
	int mtype;
	void* pt[PSIZE];
	int iparm[PSIZE];
	int maxfct, mnum, phase, error;
//	DISALLOW_COPY_AND_ASSIGN(MKLPardiso);
};

};

#endif // Interface_Pardiso_h__
