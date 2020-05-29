/*******************************************************************************
* Copyright 2005-2015 Intel Corporation All Rights Reserved.
*
* The source code,  information  and material  ("Material") contained  herein is
* owned by Intel Corporation or its  suppliers or licensors,  and  title to such
* Material remains with Intel  Corporation or its  suppliers or  licensors.  The
* Material  contains  proprietary  information  of  Intel or  its suppliers  and
* licensors.  The Material is protected by  worldwide copyright  laws and treaty
* provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
* modified, published,  uploaded, posted, transmitted,  distributed or disclosed
* in any way without Intel's prior express written permission.  No license under
* any patent,  copyright or other  intellectual property rights  in the Material
* is granted to  or  conferred  upon  you,  either   expressly,  by implication,
* inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
* property rights must be express and approved by Intel in writing.
*
* Unless otherwise agreed by Intel in writing,  you may not remove or alter this
* notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
* suppliers or licensors in any way.
*******************************************************************************/

/*
!  Content:
!   Intel(R) Math Kernel Library (MKL) interface for preconditioners, RCI ISS and
!   TR solvers routines
!******************************************************************************/

#ifndef _MKL_RCISOLVER_H_
#define _MKL_RCISOLVER_H_

#include "mkl_types.h"
#include "mkl_service.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern void dcsrilu0(MKL_INT *n, double *a,MKL_INT *ia,MKL_INT *ja, double *alu,MKL_INT *ipar, double *dpar,MKL_INT *ierr);
extern void dcsrilut(MKL_INT *n, double *a,MKL_INT *ia,MKL_INT *ja, double *alut,MKL_INT *ialut,MKL_INT *jalut,double * tol,MKL_INT *maxfil,MKL_INT *ipar, double *dpar,MKL_INT *ierr);

extern void DCSRILU0(MKL_INT *n, double *a,MKL_INT *ia,MKL_INT *ja, double *alu,MKL_INT *ipar, double *dpar,MKL_INT *ierr);
extern void DCSRILUT(MKL_INT *n, double *a,MKL_INT *ia,MKL_INT *ja, double *alut,MKL_INT *ialut,MKL_INT *jalut,double * tol,MKL_INT *maxfil,MKL_INT *ipar, double *dpar,MKL_INT *ierr);

/* PCG/PFGMRES Lower case */

extern void dcg_init(MKL_INT *n, double *x,double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void dcg_check(MKL_INT *n, double *x,double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void dcg(MKL_INT *n, double *x,double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void dcg_get(MKL_INT *n, double *x, double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp, MKL_INT *itercount);

extern void dcgmrhs_init(MKL_INT *n, double *x, MKL_INT* nRhs, double *b, MKL_INT *method, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void dcgmrhs_check(MKL_INT *n, double *x, MKL_INT* nRhs, double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void dcgmrhs(MKL_INT *n, double *x, MKL_INT* nRhs, double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void dcgmrhs_get(MKL_INT *n, double *x, MKL_INT* nRhs, double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp, MKL_INT *itercount);

extern void dfgmres_init(MKL_INT *n, double *x, double *b, MKL_INT *RCI_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void dfgmres_check(MKL_INT *n, double *x, double *b, MKL_INT *RCI_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void dfgmres(MKL_INT *n, double *x, double *b, MKL_INT *RCI_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void dfgmres_get(MKL_INT *n, double *x, double *b, MKL_INT *RCI_request, MKL_INT *ipar, double *dpar, double *tmp, MKL_INT *itercount);

/* PCG/PFGMRES Upper case */

extern void DCG_INIT(MKL_INT *n, double *x,double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void DCG_CHECK(MKL_INT *n, double *x,double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void DCG(MKL_INT *n, double *x,double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void DCG_GET(MKL_INT *n, double *x, double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp, MKL_INT *itercount);

extern void DCGMRHS_INIT(MKL_INT *n, double *x, MKL_INT* nRhs, double *b, MKL_INT *method, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void DCGMRHS_CHECK(MKL_INT *n, double *x, MKL_INT* nRhs, double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void DCGMRHS(MKL_INT *n, double *x, MKL_INT* nRhs, double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void DCGMRHS_GET(MKL_INT *n, double *x, MKL_INT* nRhs, double *b, MKL_INT *rci_request, MKL_INT *ipar, double *dpar, double *tmp, MKL_INT *itercount);

extern void DFGMRES_INIT(MKL_INT *n, double *x, double *b, MKL_INT *RCI_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void DFGMRES_CHECK(MKL_INT *n, double *x, double *b, MKL_INT *RCI_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void DFGMRES(MKL_INT *n, double *x, double *b, MKL_INT *RCI_request, MKL_INT *ipar, double *dpar, double *tmp);
extern void DFGMRES_GET(MKL_INT *n, double *x, double *b, MKL_INT *RCI_request, MKL_INT *ipar, double *dpar, double *tmp, MKL_INT *itercount);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#ifdef __cplusplus
extern "C" {
#endif

/* Return status values */
#define TR_SUCCESS        1501
#define TR_INVALID_OPTION 1502
#define TR_OUT_OF_MEMORY  1503

/* Basic data types */
typedef void* _TRNSP_HANDLE_t;
typedef void* _TRNSPBC_HANDLE_t;
typedef void* _JACOBIMATRIX_HANDLE_t;

typedef void(*USRFCND) (MKL_INT*,MKL_INT*,double*,double*);
typedef void(*USRFCNXD) (MKL_INT*,MKL_INT*,double*,double*,void*);

typedef void(*USRFCNS) (MKL_INT*,MKL_INT*,float*,float*);
typedef void(*USRFCNXS) (MKL_INT*,MKL_INT*,float*,float*,void*);

/* Function prototypes */
extern MKL_INT dtrnlsp_init     (_TRNSP_HANDLE_t*, MKL_INT*, MKL_INT*, double*, double*, MKL_INT*, MKL_INT*, double*);
extern MKL_INT dtrnlsp_check    (_TRNSP_HANDLE_t*, MKL_INT*, MKL_INT*, double*, double*, double*, MKL_INT*);
extern MKL_INT dtrnlsp_solve    (_TRNSP_HANDLE_t*, double*, double*, MKL_INT*);
extern MKL_INT dtrnlsp_get      (_TRNSP_HANDLE_t*, MKL_INT*, MKL_INT*, double*, double*);
extern MKL_INT dtrnlsp_delete   (_TRNSP_HANDLE_t*);
               
extern MKL_INT dtrnlspbc_init   (_TRNSPBC_HANDLE_t*, MKL_INT*, MKL_INT*, double*, double*, double*, double*, MKL_INT*, MKL_INT*, double*);
extern MKL_INT dtrnlspbc_check  (_TRNSPBC_HANDLE_t*, MKL_INT*, MKL_INT*, double*, double*, double*, double*, double*, MKL_INT*);
extern MKL_INT dtrnlspbc_solve  (_TRNSPBC_HANDLE_t*, double*, double*, MKL_INT*);
extern MKL_INT dtrnlspbc_get    (_TRNSPBC_HANDLE_t*, MKL_INT*, MKL_INT*, double*, double*);
extern MKL_INT dtrnlspbc_delete (_TRNSPBC_HANDLE_t*);
               
extern MKL_INT djacobi_init     (_JACOBIMATRIX_HANDLE_t*, MKL_INT*, MKL_INT*, double*, double*, double*);
extern MKL_INT djacobi_solve    (_JACOBIMATRIX_HANDLE_t*, double*, double*, MKL_INT*);
extern MKL_INT djacobi_delete   (_JACOBIMATRIX_HANDLE_t*);
extern MKL_INT djacobi          (USRFCND fcn, MKL_INT*, MKL_INT*, double*, double*, double*);
extern MKL_INT djacobix         (USRFCNXD fcn, MKL_INT*, MKL_INT*, double*, double*, double*,void*);

extern MKL_INT strnlsp_init     (_TRNSP_HANDLE_t*, MKL_INT*, MKL_INT*, float*, float*, MKL_INT*, MKL_INT*, float*);
extern MKL_INT strnlsp_check    (_TRNSP_HANDLE_t*, MKL_INT*, MKL_INT*, float*, float*, float*, MKL_INT*);
extern MKL_INT strnlsp_solve    (_TRNSP_HANDLE_t*, float*, float*, MKL_INT*);
extern MKL_INT strnlsp_get      (_TRNSP_HANDLE_t*, MKL_INT*, MKL_INT*, float*, float*);
extern MKL_INT strnlsp_delete   (_TRNSP_HANDLE_t*);
               
extern MKL_INT strnlspbc_init   (_TRNSPBC_HANDLE_t*, MKL_INT*, MKL_INT*, float*, float*, float*, float*, MKL_INT*, MKL_INT*, float*);
extern MKL_INT strnlspbc_check  (_TRNSPBC_HANDLE_t*, MKL_INT*, MKL_INT*, float*, float*, float*, float*, float*, MKL_INT*);
extern MKL_INT strnlspbc_solve  (_TRNSPBC_HANDLE_t*, float*, float*, MKL_INT*);
extern MKL_INT strnlspbc_get    (_TRNSPBC_HANDLE_t*, MKL_INT*, MKL_INT*, float*, float*);
extern MKL_INT strnlspbc_delete (_TRNSPBC_HANDLE_t*);
               
extern MKL_INT sjacobi_init     (_JACOBIMATRIX_HANDLE_t*, MKL_INT*, MKL_INT*, float*, float*, float*);
extern MKL_INT sjacobi_solve    (_JACOBIMATRIX_HANDLE_t*, float*, float*, MKL_INT*);
extern MKL_INT sjacobi_delete   (_JACOBIMATRIX_HANDLE_t*);
extern MKL_INT sjacobi          (USRFCNS fcn, MKL_INT*, MKL_INT*, float*, float*, float*);
extern MKL_INT sjacobix         (USRFCNXS fcn, MKL_INT*, MKL_INT*, float*, float*, float*,void*);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _MKL_RCISOLVER_H_ */
