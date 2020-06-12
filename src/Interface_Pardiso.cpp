#include <iostream>
#include <stdio.h>
#include "Interface_Pardiso.h"


using namespace std;

namespace UGSOL {

MKLPardiso::MKLPardiso() { 
	solution = NULL;
	phase    = 0;
	msglvl   = 0; ///< Print statistical information on screen 

    for (int i = 0; i < PSIZE; i++) {
        iparm[i] = 0;
        pt[i] = 0;
    }
}

void MKLPardiso::Clear_Pardiso() { 
	/**
	 * .. Termination and release of memory.
	 */
	if(phase!=0){
//		cout << "clear pardiso..." << endl;
        int nrhs = 1, idum;
        double ddum;
		phase = -1; ///< release internal memory.
		PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
			     &idum, &ddum, &idum, &idum, &idum, &nrhs,
			     iparm, &msglvl, &ddum, &ddum, &error);
	}

    delete[] solution;
    solution = NULL;
	phase = 0;
    msglvl= 0; 
}

MKLPardiso::~MKLPardiso() { 
	Clear_Pardiso();
}

int MKLPardiso::Init_Pardiso(int n, double* a, int* ia, int* ja, const int idx_base,
                             const int _mtype, const int _mnum, const int _cgs, const int _msglvl)
{
	mnum = _mnum;   ///< Which factorization to use. 
//	cout << "Clear Pardiso..." << endl;
	Clear_Pardiso();
//	cout << "Done" << endl;
    delete[] solution;
//	cout << "Open space for solution...";
    solution = new double[n];
//	cout << "Done" << endl;
	for(int i=0; i<n; ++i){
		solution[i] = 0.0;
	}
    int nrhs = 1, idum;
    double ddum;
	mtype = _mtype; 
//	cout << "Set up iparm...";
	iparm[0] = 1; ///< No solver default 
	iparm[1] = 2; ///< The nested dissection algorithm from the METIS package 
	iparm[2] = 1; ///< Currently is not used 
	iparm[3] = _cgs; /** If iparm[3] > 0, use CGS iterative algorithm, 
					      * taking previous LU factorization as a preconditioner. */
//	iparm[3] = 0;
	iparm[4] = 0; ///< No user fill-in reducing permutation 
	iparm[5] = 1; ///< Write solution into b 
	iparm[6] = 0; ///< Not in use 
	iparm[7] = 2; ///< Max numbers of iterative refinement steps 
	iparm[8] = 0; ///< Not in use 
	iparm[9] = 12; ///< Perturb the pivot elements with 1E-12 
	iparm[10] = 1; ///< Use nonsymmetric permutation and scaling MPS 
	iparm[11] = 0; ///< Not in use 
	iparm[12] = 1; ///< Maximum weighted matching algorithm is switched-on 
	iparm[13] = 0; ///< Output: Number of perturbed pivots 
	iparm[14] = 0; ///< Not in use 
	iparm[15] = 0; ///< Not in use 
	iparm[16] = 0; ///< Not in use 
	iparm[17] = -1; ///< Output: Number of nonzeros in the factor LU 
	iparm[18] = -1; ///< Output: Mflops for LU factorization 
	iparm[19] = 1; ///< Output: Numbers of CG Iterations 
	iparm[23] = 0; ///< Classic algorithm for factorization 

    if(idx_base==0) {
	    iparm[34] = 1;        
    } else if(idx_base==1) {
	    iparm[34] = 0;
    } else {
        printf(" Error in Init_Pardiso: illegal matrix index base.\n");
    }

	iparm[59] = 0; ///< Disallow out-of-core version 
//	cout << "Done" << endl;
	maxfct = 1;     ///< Maximum number of numerical factorizations. 
	error = 0;      ///< Initialize error flag 
//	cout << "set msg...";
	
    msglvl = _msglvl;
//	cout << msglvl << endl;
//	cout << "Done" << endl;
	/**
	 * .. Reordering and Symbolic Factorization. This step also allocates
	 * all memory that is necessary for the factorization.               
	 */
//	cout << "phase:" << phase << endl;
	//phase = 11;
	int phase_temp = 11;
//	cout << "Performing Pardiso..." ;
	//cout << maxfct <<" "<<mnum<<" "<<mtype<<" "<<phase_temp<<" "<<n<<" "<< <<endl;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase_temp,
		     &n, a, ia, ja, &idum, &nrhs,
		     iparm, &msglvl, &ddum, &ddum, &error);
//	cout << "Done..."<<endl;
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		return error;
	}
	return 0;
}

int MKLPardiso::Pardiso_Factor(int n, double* a, int* ia, int* ja){
	int tmp = iparm[3];
    int nrhs = 1, idum;
    double ddum;
	iparm[3] = 0;
	/**
	 * .. Redo-Numerical factorization, but without analyzing phase.
	 */
	phase = 22;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		     &n, a, ia, ja, &idum, &nrhs,
		     iparm, &msglvl, &ddum, &ddum, &error);

	if (error != 0) {
		printf("\nERROR during numerical factorization: %d.\n", error);
		return error;
	}
	iparm[3] = tmp;
	return 0;
}

int MKLPardiso::BackSubstitute(int n, double* a, int* ia, int* ja, double* b){
    int nrhs = 1, idum;
	phase = 33;
	PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
		     &n, a, ia, ja, &idum, &nrhs,
		     iparm, &msglvl, b, solution, &error);
	return 1;
};

int MKLPardiso::Pardiso_Solve(int n, double* a, int* ia, int* ja, double* b) 
{
    if(n>9999)
        cout<<" Direct factorization...";

    Pardiso_Factor(n, a, ia, ja);
   BackSubstitute(n, a, ia, ja, b);
    if (error != 0){
        cout<<" Pardiso error code: "<<error<<endl;
    } else {
        if(n>9999)
            cout<< " Complete." <<endl;
    }
    return error;
}

};
