/* file: mkl_vsl_functions.h */
/*******************************************************************************
* Copyright 2006-2015 Intel Corporation All Rights Reserved.
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
//++
//  User-level VSL function declarations
//--
*/

#ifndef __MKL_VSL_FUNCTIONS_H__
#define __MKL_VSL_FUNCTIONS_H__

#include "mkl_vsl_types.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*
//++
//  EXTERNAL API MACROS.
//  Used to construct VSL function declaration. Change them if you are going to
//  provide different API for VSL functions.
//--
*/

#if  !defined(_Mkl_Api)
#define _Mkl_Api(rtype,name,arg)   extern rtype name    arg;
#endif

#if  !defined(_mkl_api)
#define _mkl_api(rtype,name,arg)   extern rtype name##_ arg;
#endif

#if  !defined(_MKL_API)
#define _MKL_API(rtype,name,arg)   extern rtype name##_ arg;
#endif

/*
//++
//  VSL CONTINUOUS DISTRIBUTION GENERATOR FUNCTION DECLARATIONS.
//--
*/
/* Cauchy distribution */
_Mkl_Api(int,vdRngCauchy,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const double  , const double  ))
_MKL_API(int,VDRNGCAUCHY,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *))
_mkl_api(int,vdrngcauchy,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *))
_Mkl_Api(int,vsRngCauchy,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , float [],  const float  ,  const float   ))
_MKL_API(int,VSRNGCAUCHY,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float * ))
_mkl_api(int,vsrngcauchy,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float * ))

/* Uniform distribution */
_Mkl_Api(int,vdRngUniform,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const double  , const double  ))
_MKL_API(int,VDRNGUNIFORM,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *))
_mkl_api(int,vdrnguniform,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *))
_Mkl_Api(int,vsRngUniform,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , float [],  const float  ,  const float   ))
_MKL_API(int,VSRNGUNIFORM,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float * ))
_mkl_api(int,vsrnguniform,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float * ))

/* Gaussian distribution */
_Mkl_Api(int,vdRngGaussian,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const double  , const double  ))
_MKL_API(int,VDRNGGAUSSIAN,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *))
_mkl_api(int,vdrnggaussian,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *))
_Mkl_Api(int,vsRngGaussian,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , float [],  const float  ,  const float   ))
_MKL_API(int,VSRNGGAUSSIAN,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float * ))
_mkl_api(int,vsrnggaussian,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float * ))

/* GaussianMV distribution */
_Mkl_Api(int,vdRngGaussianMV,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const MKL_INT  ,  const MKL_INT  , const double *, const double *))
_MKL_API(int,VDRNGGAUSSIANMV,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const MKL_INT *,  const MKL_INT *, const double *, const double *))
_mkl_api(int,vdrnggaussianmv,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const MKL_INT *,  const MKL_INT *, const double *, const double *))
_Mkl_Api(int,vsRngGaussianMV,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , float [],  const MKL_INT  ,  const MKL_INT  , const float *,  const float * ))
_MKL_API(int,VSRNGGAUSSIANMV,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const MKL_INT *,  const MKL_INT *, const float *,  const float * ))
_mkl_api(int,vsrnggaussianmv,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const MKL_INT *,  const MKL_INT *, const float *,  const float * ))

/* Exponential distribution */
_Mkl_Api(int,vdRngExponential,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  ,  double [], const double  , const double  ))
_MKL_API(int,VDRNGEXPONENTIAL,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *,  double [], const double *, const double *))
_mkl_api(int,vdrngexponential,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *,  double [], const double *, const double *))
_Mkl_Api(int,vsRngExponential,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  ,  float [],  const float  ,  const float   ))
_MKL_API(int,VSRNGEXPONENTIAL,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *,  float [],  const float *,  const float * ))
_mkl_api(int,vsrngexponential,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *,  float [],  const float *,  const float * ))

/* Laplace distribution */
_Mkl_Api(int,vdRngLaplace,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const double  , const double  ))
_MKL_API(int,VDRNGLAPLACE,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *))
_mkl_api(int,vdrnglaplace,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *))
_Mkl_Api(int,vsRngLaplace,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , float [],  const float  ,  const float   ))
_MKL_API(int,VSRNGLAPLACE,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float * ))
_mkl_api(int,vsrnglaplace,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float * ))

/* Weibull distribution */
_Mkl_Api(int,vdRngWeibull,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const double  , const double  , const double  ))
_MKL_API(int,VDRNGWEIBULL,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *, const double *))
_mkl_api(int,vdrngweibull,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *, const double *))
_Mkl_Api(int,vsRngWeibull,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , float [],  const float  ,  const float  ,  const float   ))
_MKL_API(int,VSRNGWEIBULL,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float *,  const float * ))
_mkl_api(int,vsrngweibull,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float *,  const float * ))

/* Rayleigh distribution */
_Mkl_Api(int,vdRngRayleigh,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  ,  double [], const double  , const double  ))
_MKL_API(int,VDRNGRAYLEIGH,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *,  double [], const double *, const double *))
_mkl_api(int,vdrngrayleigh,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *,  double [], const double *, const double *))
_Mkl_Api(int,vsRngRayleigh,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  ,  float [],  const float  ,  const float   ))
_MKL_API(int,VSRNGRAYLEIGH,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *,  float [],  const float *,  const float * ))
_mkl_api(int,vsrngrayleigh,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *,  float [],  const float *,  const float * ))

/* Lognormal distribution */
_Mkl_Api(int,vdRngLognormal,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const double  , const double  , const double  , const double  ))
_MKL_API(int,VDRNGLOGNORMAL,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *, const double *, const double *))
_mkl_api(int,vdrnglognormal,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *, const double *, const double *))
_Mkl_Api(int,vsRngLognormal,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , float [],  const float  ,  const float  ,  const float  ,  const float   ))
_MKL_API(int,VSRNGLOGNORMAL,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float *,  const float *,  const float * ))
_mkl_api(int,vsrnglognormal,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float *,  const float *,  const float * ))

/* Gumbel distribution */
_Mkl_Api(int,vdRngGumbel,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const double  , const double  ))
_MKL_API(int,VDRNGGUMBEL,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *))
_mkl_api(int,vdrnggumbel,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *))
_Mkl_Api(int,vsRngGumbel,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , float [],  const float  ,  const float   ))
_MKL_API(int,VSRNGGUMBEL,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float * ))
_mkl_api(int,vsrnggumbel,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float * ))

/* Gamma distribution */
_Mkl_Api(int,vdRngGamma,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const double  , const double  , const double  ))
_MKL_API(int,VDRNGGAMMA,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *, const double *))
_mkl_api(int,vdrnggamma,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *, const double *))
_Mkl_Api(int,vsRngGamma,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , float [],  const float  ,  const float  ,  const float   ))
_MKL_API(int,VSRNGGAMMA,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float *,  const float * ))
_mkl_api(int,vsrnggamma,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float *,  const float * ))

/* Beta distribution */
_Mkl_Api(int,vdRngBeta,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , double [], const double  , const double  , const double  , const double  ))
_MKL_API(int,VDRNGBETA,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *, const double *, const double *))
_mkl_api(int,vdrngbeta,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, double [], const double *, const double *, const double *, const double *))
_Mkl_Api(int,vsRngBeta,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , float [],  const float  ,  const float  ,  const float  ,  const float   ))
_MKL_API(int,VSRNGBETA,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float *,  const float *,  const float * ))
_mkl_api(int,vsrngbeta,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, float [],  const float *,  const float *,  const float *,  const float * ))

/*
//++
//  VSL DISCRETE DISTRIBUTION GENERATOR FUNCTION DECLARATIONS.
//--
*/
/* Bernoulli distribution */
_Mkl_Api(int,viRngBernoulli,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , int [], const double  ))
_MKL_API(int,VIRNGBERNOULLI,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const double *))
_mkl_api(int,virngbernoulli,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const double *))

/* Uniform distribution */
_Mkl_Api(int,viRngUniform,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , int [], const int  , const int  ))
_MKL_API(int,VIRNGUNIFORM,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const int *, const int *))
_mkl_api(int,virnguniform,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const int *, const int *))

/* UniformBits distribution */
_Mkl_Api(int,viRngUniformBits,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , unsigned int []))
_MKL_API(int,VIRNGUNIFORMBITS,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, unsigned int []))
_mkl_api(int,virnguniformbits,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, unsigned int []))

/* UniformBits32 distribution */
_Mkl_Api(int,viRngUniformBits32,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , unsigned int []))
_MKL_API(int,VIRNGUNIFORMBITS32,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, unsigned int []))
_mkl_api(int,virnguniformbits32,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, unsigned int []))

/* UniformBits64 distribution */
_Mkl_Api(int,viRngUniformBits64,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , unsigned MKL_INT64 []))
_MKL_API(int,VIRNGUNIFORMBITS64,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, unsigned MKL_INT64 []))
_mkl_api(int,virnguniformbits64,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, unsigned MKL_INT64 []))

/* Geometric distribution */
_Mkl_Api(int,viRngGeometric,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , int [], const double  ))
_MKL_API(int,VIRNGGEOMETRIC,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const double *))
_mkl_api(int,virnggeometric,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const double *))

/* Binomial distribution */
_Mkl_Api(int,viRngBinomial,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , int [], const int  , const double  ))
_MKL_API(int,VIRNGBINOMIAL,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const int *, const double *))
_mkl_api(int,virngbinomial,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const int *, const double *))

/* Hypergeometric distribution */
_Mkl_Api(int,viRngHypergeometric,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , int [], const int  , const int  , const int  ))
_MKL_API(int,VIRNGHYPERGEOMETRIC,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const int *, const int *, const int *))
_mkl_api(int,virnghypergeometric,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const int *, const int *, const int *))

/* Negbinomial distribution */
_Mkl_Api(int,viRngNegbinomial,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , int [], const double  , const double  ))
_Mkl_Api(int,viRngNegBinomial,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , int [], const double  , const double  ))
_MKL_API(int,VIRNGNEGBINOMIAL,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const double *, const double *))
_mkl_api(int,virngnegbinomial,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const double *, const double *))

/* Poisson distribution */
_Mkl_Api(int,viRngPoisson,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , int [], const double  ))
_MKL_API(int,VIRNGPOISSON,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const double *))
_mkl_api(int,virngpoisson,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const double *))

/* PoissonV distribution */
_Mkl_Api(int,viRngPoissonV,(const MKL_INT  , VSLStreamStatePtr  , const MKL_INT  , int [], const double []))
_MKL_API(int,VIRNGPOISSONV,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const double []))
_mkl_api(int,virngpoissonv,(const MKL_INT *, VSLStreamStatePtr *, const MKL_INT *, int [], const double []))


/*
//++
//  VSL SERVICE FUNCTION DECLARATIONS.
//--
*/
/* NewStream - stream creation/initialization */
_Mkl_Api(int,vslNewStream,(VSLStreamStatePtr* , const MKL_INT  , const MKL_UINT  ))
_mkl_api(int,vslnewstream,(VSLStreamStatePtr* , const MKL_INT *, const MKL_UINT *))
_MKL_API(int,VSLNEWSTREAM,(VSLStreamStatePtr* , const MKL_INT *, const MKL_UINT *))

/* NewStreamEx - advanced stream creation/initialization */
_Mkl_Api(int,vslNewStreamEx,(VSLStreamStatePtr* , const MKL_INT  , const MKL_INT  , const unsigned int[]))
_mkl_api(int,vslnewstreamex,(VSLStreamStatePtr* , const MKL_INT *, const MKL_INT *, const unsigned int[]))
_MKL_API(int,VSLNEWSTREAMEX,(VSLStreamStatePtr* , const MKL_INT *, const MKL_INT *, const unsigned int[]))

_Mkl_Api(int,vsliNewAbstractStream,(VSLStreamStatePtr* , const MKL_INT  , const unsigned int[], const iUpdateFuncPtr))
_mkl_api(int,vslinewabstractstream,(VSLStreamStatePtr* , const MKL_INT *, const unsigned int[], const iUpdateFuncPtr))
_MKL_API(int,VSLINEWABSTRACTSTREAM,(VSLStreamStatePtr* , const MKL_INT *, const unsigned int[], const iUpdateFuncPtr))

_Mkl_Api(int,vsldNewAbstractStream,(VSLStreamStatePtr* , const MKL_INT  , const double[], const double  , const double  , const dUpdateFuncPtr))
_mkl_api(int,vsldnewabstractstream,(VSLStreamStatePtr* , const MKL_INT *, const double[], const double *, const double *, const dUpdateFuncPtr))
_MKL_API(int,VSLDNEWABSTRACTSTREAM,(VSLStreamStatePtr* , const MKL_INT *, const double[], const double *, const double *, const dUpdateFuncPtr))

_Mkl_Api(int,vslsNewAbstractStream,(VSLStreamStatePtr* , const MKL_INT  , const float[], const float  , const float  , const sUpdateFuncPtr))
_mkl_api(int,vslsnewabstractstream,(VSLStreamStatePtr* , const MKL_INT *, const float[], const float *, const float *, const sUpdateFuncPtr))
_MKL_API(int,VSLSNEWABSTRACTSTREAM,(VSLStreamStatePtr* , const MKL_INT *, const float[], const float *, const float *, const sUpdateFuncPtr))

/* DeleteStream - delete stream */
_Mkl_Api(int,vslDeleteStream,(VSLStreamStatePtr*))
_mkl_api(int,vsldeletestream,(VSLStreamStatePtr*))
_MKL_API(int,VSLDELETESTREAM,(VSLStreamStatePtr*))

/* CopyStream - copy all stream information */
_Mkl_Api(int,vslCopyStream,(VSLStreamStatePtr*, const VSLStreamStatePtr))
_mkl_api(int,vslcopystream,(VSLStreamStatePtr*, const VSLStreamStatePtr))
_MKL_API(int,VSLCOPYSTREAM,(VSLStreamStatePtr*, const VSLStreamStatePtr))

/* CopyStreamState - copy stream state only */
_Mkl_Api(int,vslCopyStreamState,(VSLStreamStatePtr  , const VSLStreamStatePtr  ))
_mkl_api(int,vslcopystreamstate,(VSLStreamStatePtr *, const VSLStreamStatePtr *))
_MKL_API(int,VSLCOPYSTREAMSTATE,(VSLStreamStatePtr *, const VSLStreamStatePtr *))

/* LeapfrogStream - leapfrog method */
_Mkl_Api(int,vslLeapfrogStream,(VSLStreamStatePtr  , const MKL_INT  , const MKL_INT  ))
_mkl_api(int,vslleapfrogstream,(VSLStreamStatePtr *, const MKL_INT *, const MKL_INT *))
_MKL_API(int,VSLLEAPFROGSTREAM,(VSLStreamStatePtr *, const MKL_INT *, const MKL_INT *))

/* SkipAheadStream - skip-ahead method */
_Mkl_Api(int,vslSkipAheadStream,(VSLStreamStatePtr  , const long long int  ))
_mkl_api(int,vslskipaheadstream,(VSLStreamStatePtr *, const long long int *))
_MKL_API(int,VSLSKIPAHEADSTREAM,(VSLStreamStatePtr *, const long long int *))

/* GetStreamStateBrng - get BRNG associated with given stream */
_Mkl_Api(int,vslGetStreamStateBrng,(const VSLStreamStatePtr  ))
_mkl_api(int,vslgetstreamstatebrng,(const VSLStreamStatePtr *))
_MKL_API(int,VSLGETSTREAMSTATEBRNG,(const VSLStreamStatePtr *))

/* GetNumRegBrngs - get number of registered BRNGs */
_Mkl_Api(int,vslGetNumRegBrngs,(void))
_mkl_api(int,vslgetnumregbrngs,(void))
_MKL_API(int,VSLGETNUMREGBRNGS,(void))

/* RegisterBrng - register new BRNG */
_Mkl_Api(int,vslRegisterBrng,(const VSLBRngProperties* ))
_mkl_api(int,vslregisterbrng,(const VSLBRngProperties* ))
_MKL_API(int,VSLREGISTERBRNG,(const VSLBRngProperties* ))

/* GetBrngProperties - get BRNG properties */
_Mkl_Api(int,vslGetBrngProperties,(const int  , VSLBRngProperties* ))
_mkl_api(int,vslgetbrngproperties,(const int *, VSLBRngProperties* ))
_MKL_API(int,VSLGETBRNGPROPERTIES,(const int *, VSLBRngProperties* ))

/* SaveStreamF - save random stream descriptive data to file */
_Mkl_Api(int,vslSaveStreamF,(const VSLStreamStatePtr  , const char*             ))
_mkl_api(int,vslsavestreamf,(const VSLStreamStatePtr *, const char* , const int ))
_MKL_API(int,VSLSAVESTREAMF,(const VSLStreamStatePtr *, const char* , const int ))

/* LoadStreamF - load random stream descriptive data from file */
_Mkl_Api(int,vslLoadStreamF,(VSLStreamStatePtr *, const char*             ))
_mkl_api(int,vslloadstreamf,(VSLStreamStatePtr *, const char* , const int ))
_MKL_API(int,VSLLOADSTREAMF,(VSLStreamStatePtr *, const char* , const int ))

/* SaveStreamM - save random stream descriptive data to memory */
_Mkl_Api(int,vslSaveStreamM,(const VSLStreamStatePtr  , char* ))
_mkl_api(int,vslsavestreamm,(const VSLStreamStatePtr *, char* ))
_MKL_API(int,VSLSAVESTREAMM,(const VSLStreamStatePtr *, char* ))

/* LoadStreamM - load random stream descriptive data from memory */
_Mkl_Api(int,vslLoadStreamM,(VSLStreamStatePtr *, const char* ))
_mkl_api(int,vslloadstreamm,(VSLStreamStatePtr *, const char* ))
_MKL_API(int,VSLLOADSTREAMM,(VSLStreamStatePtr *, const char* ))

/* GetStreamSize - get size of random stream descriptive data */
_Mkl_Api(int,vslGetStreamSize,(const VSLStreamStatePtr))
_mkl_api(int,vslgetstreamsize,(const VSLStreamStatePtr))
_MKL_API(int,VSLGETSTREAMSIZE,(const VSLStreamStatePtr))

/*
//++
//  VSL CONVOLUTION AND CORRELATION FUNCTION DECLARATIONS.
//--
*/

_Mkl_Api(int,vsldConvNewTask,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_mkl_api(int,vsldconvnewtask,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_MKL_API(int,VSLDCONVNEWTASK,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))

_Mkl_Api(int,vslsConvNewTask,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_mkl_api(int,vslsconvnewtask,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_MKL_API(int,VSLSCONVNEWTASK,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))

_Mkl_Api(int,vslzConvNewTask,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_mkl_api(int,vslzconvnewtask,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_MKL_API(int,VSLZCONVNEWTASK,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))

_Mkl_Api(int,vslcConvNewTask,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_mkl_api(int,vslcconvnewtask,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_MKL_API(int,VSLCCONVNEWTASK,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))

_Mkl_Api(int,vsldCorrNewTask,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_mkl_api(int,vsldcorrnewtask,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_MKL_API(int,VSLDCORRNEWTASK,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))

_Mkl_Api(int,vslsCorrNewTask,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_mkl_api(int,vslscorrnewtask,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_MKL_API(int,VSLSCORRNEWTASK,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))

_Mkl_Api(int,vslzCorrNewTask,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_mkl_api(int,vslzcorrnewtask,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_MKL_API(int,VSLZCORRNEWTASK,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))

_Mkl_Api(int,vslcCorrNewTask,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_mkl_api(int,vslccorrnewtask,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))
_MKL_API(int,VSLCCORRNEWTASK,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT []))


_Mkl_Api(int,vsldConvNewTask1D,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT ,  const MKL_INT  ))
_mkl_api(int,vsldconvnewtask1d,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))
_MKL_API(int,VSLDCONVNEWTASK1D,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))

_Mkl_Api(int,vslsConvNewTask1D,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT ,  const MKL_INT  ))
_mkl_api(int,vslsconvnewtask1d,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))
_MKL_API(int,VSLSCONVNEWTASK1D,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))

_Mkl_Api(int,vslzConvNewTask1D,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT ,  const MKL_INT  ))
_mkl_api(int,vslzconvnewtask1d,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))
_MKL_API(int,VSLZCONVNEWTASK1D,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))

_Mkl_Api(int,vslcConvNewTask1D,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT ,  const MKL_INT  ))
_mkl_api(int,vslcconvnewtask1d,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))
_MKL_API(int,VSLCCONVNEWTASK1D,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))

_Mkl_Api(int,vsldCorrNewTask1D,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT ,  const MKL_INT  ))
_mkl_api(int,vsldcorrnewtask1d,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))
_MKL_API(int,VSLDCORRNEWTASK1D,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))

_Mkl_Api(int,vslsCorrNewTask1D,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT ,  const MKL_INT  ))
_mkl_api(int,vslscorrnewtask1d,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))
_MKL_API(int,VSLSCORRNEWTASK1D,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))

_Mkl_Api(int,vslzCorrNewTask1D,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT ,  const MKL_INT  ))
_mkl_api(int,vslzcorrnewtask1d,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))
_MKL_API(int,VSLZCORRNEWTASK1D,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))

_Mkl_Api(int,vslcCorrNewTask1D,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT ,  const MKL_INT  ))
_mkl_api(int,vslccorrnewtask1d,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))
_MKL_API(int,VSLCCORRNEWTASK1D,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* ))


_Mkl_Api(int,vsldConvNewTaskX,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT [], const double [], const MKL_INT []))
_mkl_api(int,vsldconvnewtaskx,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const double [], const MKL_INT []))
_MKL_API(int,VSLDCONVNEWTASKX,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const double [], const MKL_INT []))

_Mkl_Api(int,vslsConvNewTaskX,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT [], const float [],  const MKL_INT []))
_mkl_api(int,vslsconvnewtaskx,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const float [],  const MKL_INT []))
_MKL_API(int,VSLSCONVNEWTASKX,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const float [],  const MKL_INT []))

_Mkl_Api(int,vslzConvNewTaskX,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT []))
_mkl_api(int,vslzconvnewtaskx,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT []))
_MKL_API(int,VSLZCONVNEWTASKX,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT []))

_Mkl_Api(int,vslcConvNewTaskX,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex8 [],  const MKL_INT []))
_mkl_api(int,vslcconvnewtaskx,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex8 [],  const MKL_INT []))
_MKL_API(int,VSLCCONVNEWTASKX,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex8 [],  const MKL_INT []))

_Mkl_Api(int,vsldCorrNewTaskX,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT [], const double [], const MKL_INT []))
_mkl_api(int,vsldcorrnewtaskx,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const double [], const MKL_INT []))
_MKL_API(int,VSLDCORRNEWTASKX,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const double [], const MKL_INT []))

_Mkl_Api(int,vslsCorrNewTaskX,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT [], const float [],  const MKL_INT []))
_mkl_api(int,vslscorrnewtaskx,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const float [],  const MKL_INT []))
_MKL_API(int,VSLSCORRNEWTASKX,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const float [],  const MKL_INT []))

_Mkl_Api(int,vslzCorrNewTaskX,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT []))
_mkl_api(int,vslzcorrnewtaskx,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT []))
_MKL_API(int,VSLZCORRNEWTASKX,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT []))

_Mkl_Api(int,vslcCorrNewTaskX,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex8 [],  const MKL_INT []))
_mkl_api(int,vslccorrnewtaskx,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex8 [],  const MKL_INT []))
_MKL_API(int,VSLCCORRNEWTASKX,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT [], const MKL_INT [], const MKL_INT [], const MKL_Complex8 [],  const MKL_INT []))


_Mkl_Api(int,vsldConvNewTaskX1D,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_INT  , const double [], const MKL_INT  ))
_mkl_api(int,vsldconvnewtaskx1d,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const double [], const MKL_INT* ))
_MKL_API(int,VSLDCONVNEWTASKX1D,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const double [], const MKL_INT* ))

_Mkl_Api(int,vslsConvNewTaskX1D,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_INT  , const float [],  const MKL_INT  ))
_mkl_api(int,vslsconvnewtaskx1d,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const float [],  const MKL_INT* ))
_MKL_API(int,VSLSCONVNEWTASKX1D,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const float [],  const MKL_INT* ))

_Mkl_Api(int,vslzConvNewTaskX1D,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_Complex16 [], const MKL_INT  ))
_mkl_api(int,vslzconvnewtaskx1d,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_Complex16 [], const MKL_INT* ))
_MKL_API(int,VSLZCONVNEWTASKX1D,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_Complex16 [], const MKL_INT* ))

_Mkl_Api(int,vslcConvNewTaskX1D,(VSLConvTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_Complex8 [],  const MKL_INT  ))
_mkl_api(int,vslcconvnewtaskx1d,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_Complex8 [],  const MKL_INT* ))
_MKL_API(int,VSLCCONVNEWTASKX1D,(VSLConvTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_Complex8 [],  const MKL_INT* ))

_Mkl_Api(int,vsldCorrNewTaskX1D,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_INT  , const double [], const MKL_INT  ))
_mkl_api(int,vsldcorrnewtaskx1d,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const double [], const MKL_INT* ))
_MKL_API(int,VSLDCORRNEWTASKX1D,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const double [], const MKL_INT* ))

_Mkl_Api(int,vslsCorrNewTaskX1D,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_INT  , const float [],  const MKL_INT  ))
_mkl_api(int,vslscorrnewtaskx1d,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const float [],  const MKL_INT* ))
_MKL_API(int,VSLSCORRNEWTASKX1D,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const float [],  const MKL_INT* ))

_Mkl_Api(int,vslzCorrNewTaskX1D,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_Complex16 [], const MKL_INT  ))
_mkl_api(int,vslzcorrnewtaskx1d,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_Complex16 [], const MKL_INT* ))
_MKL_API(int,VSLZCORRNEWTASKX1D,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_Complex16 [], const MKL_INT* ))

_Mkl_Api(int,vslcCorrNewTaskX1D,(VSLCorrTaskPtr* , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_INT  , const MKL_Complex8 [],  const MKL_INT  ))
_mkl_api(int,vslccorrnewtaskx1d,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_Complex8 [],  const MKL_INT* ))
_MKL_API(int,VSLCCORRNEWTASKX1D,(VSLCorrTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const MKL_Complex8 [],  const MKL_INT* ))


_Mkl_Api(int,vslConvDeleteTask,(VSLConvTaskPtr* ))
_mkl_api(int,vslconvdeletetask,(VSLConvTaskPtr* ))
_MKL_API(int,VSLCONVDeleteTask,(VSLConvTaskPtr* ))

_Mkl_Api(int,vslCorrDeleteTask,(VSLCorrTaskPtr* ))
_mkl_api(int,vslcorrdeletetask,(VSLCorrTaskPtr* ))
_MKL_API(int,VSLCORRDeleteTask,(VSLCorrTaskPtr* ))


_Mkl_Api(int,vslConvCopyTask,(VSLConvTaskPtr* , const VSLConvTaskPtr  ))
_mkl_api(int,vslconvcopytask,(VSLConvTaskPtr* , const VSLConvTaskPtr* ))
_MKL_API(int,VSLCONVCopyTask,(VSLConvTaskPtr* , const VSLConvTaskPtr* ))

_Mkl_Api(int,vslCorrCopyTask,(VSLCorrTaskPtr* , const VSLCorrTaskPtr  ))
_mkl_api(int,vslcorrcopytask,(VSLCorrTaskPtr* , const VSLCorrTaskPtr* ))
_MKL_API(int,VSLCORRCopyTask,(VSLCorrTaskPtr* , const VSLCorrTaskPtr* ))


_Mkl_Api(int,vslConvSetMode,(VSLConvTaskPtr  , const MKL_INT  ))
_mkl_api(int,vslconvsetmode,(VSLConvTaskPtr* , const MKL_INT* ))
_MKL_API(int,VSLCONVSETMODE,(VSLConvTaskPtr* , const MKL_INT* ))

_Mkl_Api(int,vslCorrSetMode,(VSLCorrTaskPtr  , const MKL_INT  ))
_mkl_api(int,vslcorrsetmode,(VSLCorrTaskPtr* , const MKL_INT* ))
_MKL_API(int,VSLCORRSETMODE,(VSLCorrTaskPtr* , const MKL_INT* ))


_Mkl_Api(int,vslConvSetInternalPrecision,(VSLConvTaskPtr  , const MKL_INT  ))
_mkl_api(int,vslconvsetinternalprecision,(VSLConvTaskPtr* , const MKL_INT* ))
_MKL_API(int,VSLCONVSETINTERNALPRECISION,(VSLConvTaskPtr* , const MKL_INT* ))

_Mkl_Api(int,vslCorrSetInternalPrecision,(VSLCorrTaskPtr  , const MKL_INT  ))
_mkl_api(int,vslcorrsetinternalprecision,(VSLCorrTaskPtr* , const MKL_INT* ))
_MKL_API(int,VSLCORRSETINTERNALPRECISION,(VSLCorrTaskPtr* , const MKL_INT* ))


_Mkl_Api(int,vslConvSetStart,(VSLConvTaskPtr  , const MKL_INT []))
_mkl_api(int,vslconvsetstart,(VSLConvTaskPtr* , const MKL_INT []))
_MKL_API(int,VSLCONVSETSTART,(VSLConvTaskPtr* , const MKL_INT []))

_Mkl_Api(int,vslCorrSetStart,(VSLCorrTaskPtr  , const MKL_INT []))
_mkl_api(int,vslcorrsetstart,(VSLCorrTaskPtr* , const MKL_INT []))
_MKL_API(int,VSLCORRSETSTART,(VSLCorrTaskPtr* , const MKL_INT []))


_Mkl_Api(int,vslConvSetDecimation,(VSLConvTaskPtr  , const MKL_INT []))
_mkl_api(int,vslconvsetdecimation,(VSLConvTaskPtr* , const MKL_INT []))
_MKL_API(int,VSLCONVSETDECIMATION,(VSLConvTaskPtr* , const MKL_INT []))

_Mkl_Api(int,vslCorrSetDecimation,(VSLCorrTaskPtr  , const MKL_INT []))
_mkl_api(int,vslcorrsetdecimation,(VSLCorrTaskPtr* , const MKL_INT []))
_MKL_API(int,VSLCORRSETDECIMATION,(VSLCorrTaskPtr* , const MKL_INT []))


_Mkl_Api(int,vsldConvExec,(VSLConvTaskPtr  , const double [], const MKL_INT [], const double [], const MKL_INT [], double [], const MKL_INT []))
_mkl_api(int,vsldconvexec,(VSLConvTaskPtr* , const double [], const MKL_INT [], const double [], const MKL_INT [], double [], const MKL_INT []))
_MKL_API(int,VSLDCONVEXEC,(VSLConvTaskPtr* , const double [], const MKL_INT [], const double [], const MKL_INT [], double [], const MKL_INT []))

_Mkl_Api(int,vslsConvExec,(VSLConvTaskPtr  , const float [],  const MKL_INT [], const float [],  const MKL_INT [], float [],  const MKL_INT []))
_mkl_api(int,vslsconvexec,(VSLConvTaskPtr* , const float [],  const MKL_INT [], const float [],  const MKL_INT [], float [],  const MKL_INT []))
_MKL_API(int,VSLSCONVEXEC,(VSLConvTaskPtr* , const float [],  const MKL_INT [], const float [],  const MKL_INT [], float [],  const MKL_INT []))

_Mkl_Api(int,vslzConvExec,(VSLConvTaskPtr  , const MKL_Complex16 [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))
_mkl_api(int,vslzconvexec,(VSLConvTaskPtr* , const MKL_Complex16 [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))
_MKL_API(int,VSLZCONVEXEC,(VSLConvTaskPtr* , const MKL_Complex16 [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))

_Mkl_Api(int,vslcConvExec,(VSLConvTaskPtr  , const MKL_Complex8 [],  const MKL_INT [], const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))
_mkl_api(int,vslcconvexec,(VSLConvTaskPtr* , const MKL_Complex8 [],  const MKL_INT [], const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))
_MKL_API(int,VSLCCONVEXEC,(VSLConvTaskPtr* , const MKL_Complex8 [],  const MKL_INT [], const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))

_Mkl_Api(int,vsldCorrExec,(VSLCorrTaskPtr  , const double [], const MKL_INT [], const double [], const MKL_INT [], double [], const MKL_INT []))
_mkl_api(int,vsldcorrexec,(VSLCorrTaskPtr* , const double [], const MKL_INT [], const double [], const MKL_INT [], double [], const MKL_INT []))
_MKL_API(int,VSLDCORREXEC,(VSLCorrTaskPtr* , const double [], const MKL_INT [], const double [], const MKL_INT [], double [], const MKL_INT []))

_Mkl_Api(int,vslsCorrExec,(VSLCorrTaskPtr  , const float [],  const MKL_INT [], const float [],  const MKL_INT [], float [],  const MKL_INT []))
_mkl_api(int,vslscorrexec,(VSLCorrTaskPtr* , const float [],  const MKL_INT [], const float [],  const MKL_INT [], float [],  const MKL_INT []))
_MKL_API(int,VSLSCORREXEC,(VSLCorrTaskPtr* , const float [],  const MKL_INT [], const float [],  const MKL_INT [], float [],  const MKL_INT []))

_Mkl_Api(int,vslzCorrExec,(VSLCorrTaskPtr  , const MKL_Complex16 [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))
_mkl_api(int,vslzcorrexec,(VSLCorrTaskPtr* , const MKL_Complex16 [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))
_MKL_API(int,VSLZCORREXEC,(VSLCorrTaskPtr* , const MKL_Complex16 [], const MKL_INT [], const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))

_Mkl_Api(int,vslcCorrExec,(VSLCorrTaskPtr  , const MKL_Complex8 [],  const MKL_INT [], const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))
_mkl_api(int,vslccorrexec,(VSLCorrTaskPtr* , const MKL_Complex8 [],  const MKL_INT [], const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))
_MKL_API(int,VSLCCORREXEC,(VSLCorrTaskPtr* , const MKL_Complex8 [],  const MKL_INT [], const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))


_Mkl_Api(int,vsldConvExec1D,(VSLConvTaskPtr  , const double [], const MKL_INT  , const double [], const MKL_INT  , double [], const MKL_INT  ))
_mkl_api(int,vsldconvexec1d,(VSLConvTaskPtr* , const double [], const MKL_INT* , const double [], const MKL_INT* , double [], const MKL_INT* ))
_MKL_API(int,VSLDCONVEXEC1D,(VSLConvTaskPtr* , const double [], const MKL_INT* , const double [], const MKL_INT* , double [], const MKL_INT* ))

_Mkl_Api(int,vslsConvExec1D,(VSLConvTaskPtr  , const float [],  const MKL_INT  , const float [],  const MKL_INT  , float [],  const MKL_INT  ))
_mkl_api(int,vslsconvexec1d,(VSLConvTaskPtr* , const float [],  const MKL_INT* , const float [],  const MKL_INT* , float [],  const MKL_INT* ))
_MKL_API(int,VSLSCONVEXEC1D,(VSLConvTaskPtr* , const float [],  const MKL_INT* , const float [],  const MKL_INT* , float [],  const MKL_INT* ))

_Mkl_Api(int,vslzConvExec1D,(VSLConvTaskPtr  , const MKL_Complex16 [], const MKL_INT  , const MKL_Complex16 [], const MKL_INT  , MKL_Complex16 [], const MKL_INT  ))
_mkl_api(int,vslzconvexec1d,(VSLConvTaskPtr* , const MKL_Complex16 [], const MKL_INT* , const MKL_Complex16 [], const MKL_INT* , MKL_Complex16 [], const MKL_INT* ))
_MKL_API(int,VSLZCONVEXEC1D,(VSLConvTaskPtr* , const MKL_Complex16 [], const MKL_INT* , const MKL_Complex16 [], const MKL_INT* , MKL_Complex16 [], const MKL_INT* ))

_Mkl_Api(int,vslcConvExec1D,(VSLConvTaskPtr  , const MKL_Complex8 [],  const MKL_INT  , const MKL_Complex8 [],  const MKL_INT  , MKL_Complex8 [],  const MKL_INT  ))
_mkl_api(int,vslcconvexec1d,(VSLConvTaskPtr* , const MKL_Complex8 [],  const MKL_INT* , const MKL_Complex8 [],  const MKL_INT* , MKL_Complex8 [],  const MKL_INT* ))
_MKL_API(int,VSLCCONVEXEC1D,(VSLConvTaskPtr* , const MKL_Complex8 [],  const MKL_INT* , const MKL_Complex8 [],  const MKL_INT* , MKL_Complex8 [],  const MKL_INT* ))

_Mkl_Api(int,vsldCorrExec1D,(VSLCorrTaskPtr  , const double [], const MKL_INT  , const double [], const MKL_INT  , double [], const MKL_INT  ))
_mkl_api(int,vsldcorrexec1d,(VSLCorrTaskPtr* , const double [], const MKL_INT* , const double [], const MKL_INT* , double [], const MKL_INT* ))
_MKL_API(int,VSLDCORREXEC1D,(VSLCorrTaskPtr* , const double [], const MKL_INT* , const double [], const MKL_INT* , double [], const MKL_INT* ))

_Mkl_Api(int,vslsCorrExec1D,(VSLCorrTaskPtr  , const float [],  const MKL_INT  , const float [],  const MKL_INT  , float [],  const MKL_INT  ))
_mkl_api(int,vslscorrexec1d,(VSLCorrTaskPtr* , const float [],  const MKL_INT* , const float [],  const MKL_INT* , float [],  const MKL_INT* ))
_MKL_API(int,VSLSCORREXEC1D,(VSLCorrTaskPtr* , const float [],  const MKL_INT* , const float [],  const MKL_INT* , float [],  const MKL_INT* ))

_Mkl_Api(int,vslzCorrExec1D,(VSLCorrTaskPtr  , const MKL_Complex16 [], const MKL_INT  , const MKL_Complex16 [], const MKL_INT  , MKL_Complex16 [], const MKL_INT  ))
_mkl_api(int,vslzcorrexec1d,(VSLCorrTaskPtr* , const MKL_Complex16 [], const MKL_INT* , const MKL_Complex16 [], const MKL_INT* , MKL_Complex16 [], const MKL_INT* ))
_MKL_API(int,VSLZCORREXEC1D,(VSLCorrTaskPtr* , const MKL_Complex16 [], const MKL_INT* , const MKL_Complex16 [], const MKL_INT* , MKL_Complex16 [], const MKL_INT* ))

_Mkl_Api(int,vslcCorrExec1D,(VSLCorrTaskPtr  , const MKL_Complex8 [],  const MKL_INT  , const MKL_Complex8 [],  const MKL_INT  , MKL_Complex8 [],  const MKL_INT  ))
_mkl_api(int,vslccorrexec1d,(VSLCorrTaskPtr* , const MKL_Complex8 [],  const MKL_INT* , const MKL_Complex8 [],  const MKL_INT* , MKL_Complex8 [],  const MKL_INT* ))
_MKL_API(int,VSLCCORREXEC1D,(VSLCorrTaskPtr* , const MKL_Complex8 [],  const MKL_INT* , const MKL_Complex8 [],  const MKL_INT* , MKL_Complex8 [],  const MKL_INT* ))


_Mkl_Api(int,vsldConvExecX,(VSLConvTaskPtr  , const double [], const MKL_INT [], double [], const MKL_INT []))
_mkl_api(int,vsldconvexecx,(VSLConvTaskPtr* , const double [], const MKL_INT [], double [], const MKL_INT []))
_MKL_API(int,VSLDCONVEXECX,(VSLConvTaskPtr* , const double [], const MKL_INT [], double [], const MKL_INT []))

_Mkl_Api(int,vslsConvExecX,(VSLConvTaskPtr  , const float [],  const MKL_INT [], float [],  const MKL_INT []))
_mkl_api(int,vslsconvexecx,(VSLConvTaskPtr* , const float [],  const MKL_INT [], float [],  const MKL_INT []))
_MKL_API(int,VSLSCONVEXECX,(VSLConvTaskPtr* , const float [],  const MKL_INT [], float [],  const MKL_INT []))

_Mkl_Api(int,vslzConvExecX,(VSLConvTaskPtr  , const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))
_mkl_api(int,vslzconvexecx,(VSLConvTaskPtr* , const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))
_MKL_API(int,VSLZCONVEXECX,(VSLConvTaskPtr* , const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))

_Mkl_Api(int,vslcConvExecX,(VSLConvTaskPtr  , const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))
_mkl_api(int,vslcconvexecx,(VSLConvTaskPtr* , const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))
_MKL_API(int,VSLCCONVEXECX,(VSLConvTaskPtr* , const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))

_Mkl_Api(int,vsldCorrExecX,(VSLCorrTaskPtr  , const double [], const MKL_INT [], double [], const MKL_INT []))
_mkl_api(int,vsldcorrexecx,(VSLCorrTaskPtr* , const double [], const MKL_INT [], double [], const MKL_INT []))
_MKL_API(int,VSLDCORREXECX,(VSLCorrTaskPtr* , const double [], const MKL_INT [], double [], const MKL_INT []))

_Mkl_Api(int,vslsCorrExecX,(VSLCorrTaskPtr  , const float [],  const MKL_INT [], float [],  const MKL_INT []))
_mkl_api(int,vslscorrexecx,(VSLCorrTaskPtr* , const float [],  const MKL_INT [], float [],  const MKL_INT []))
_MKL_API(int,VSLSCORREXECX,(VSLCorrTaskPtr* , const float [],  const MKL_INT [], float [],  const MKL_INT []))

_Mkl_Api(int,vslzCorrExecX,(VSLCorrTaskPtr  , const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))
_mkl_api(int,vslzcorrexecx,(VSLCorrTaskPtr* , const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))
_MKL_API(int,VSLZCORREXECX,(VSLCorrTaskPtr* , const MKL_Complex16 [], const MKL_INT [], MKL_Complex16 [], const MKL_INT []))

_Mkl_Api(int,vslcCorrExecX,(VSLCorrTaskPtr  , const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))
_mkl_api(int,vslccorrexecx,(VSLCorrTaskPtr* , const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))
_MKL_API(int,VSLCCORREXECX,(VSLCorrTaskPtr* , const MKL_Complex8 [],  const MKL_INT [], MKL_Complex8 [],  const MKL_INT []))


_Mkl_Api(int,vsldConvExecX1D,(VSLConvTaskPtr  , const double [], const MKL_INT  , double [], const MKL_INT  ))
_mkl_api(int,vsldconvexecx1d,(VSLConvTaskPtr* , const double [], const MKL_INT* , double [], const MKL_INT* ))
_MKL_API(int,VSLDCONVEXECX1D,(VSLConvTaskPtr* , const double [], const MKL_INT* , double [], const MKL_INT* ))

_Mkl_Api(int,vslsConvExecX1D,(VSLConvTaskPtr  , const float [],  const MKL_INT  , float [],  const MKL_INT  ))
_mkl_api(int,vslsconvexecx1d,(VSLConvTaskPtr* , const float [],  const MKL_INT* , float [],  const MKL_INT* ))
_MKL_API(int,VSLSCONVEXECX1D,(VSLConvTaskPtr* , const float [],  const MKL_INT* , float [],  const MKL_INT* ))

_Mkl_Api(int,vslzConvExecX1D,(VSLConvTaskPtr  , const MKL_Complex16 [], const MKL_INT  , MKL_Complex16 [], const MKL_INT  ))
_mkl_api(int,vslzconvexecx1d,(VSLConvTaskPtr* , const MKL_Complex16 [], const MKL_INT* , MKL_Complex16 [], const MKL_INT* ))
_MKL_API(int,VSLZCONVEXECX1D,(VSLConvTaskPtr* , const MKL_Complex16 [], const MKL_INT* , MKL_Complex16 [], const MKL_INT* ))

_Mkl_Api(int,vslcConvExecX1D,(VSLConvTaskPtr  , const MKL_Complex8 [],  const MKL_INT  , MKL_Complex8 [],  const MKL_INT  ))
_mkl_api(int,vslcconvexecx1d,(VSLConvTaskPtr* , const MKL_Complex8 [],  const MKL_INT* , MKL_Complex8 [],  const MKL_INT* ))
_MKL_API(int,VSLCCONVEXECX1D,(VSLConvTaskPtr* , const MKL_Complex8 [],  const MKL_INT* , MKL_Complex8 [],  const MKL_INT* ))

_Mkl_Api(int,vsldCorrExecX1D,(VSLCorrTaskPtr  , const double [], const MKL_INT  , double [], const MKL_INT  ))
_mkl_api(int,vsldcorrexecx1d,(VSLCorrTaskPtr* , const double [], const MKL_INT* , double [], const MKL_INT* ))
_MKL_API(int,VSLDCORREXECX1D,(VSLCorrTaskPtr* , const double [], const MKL_INT* , double [], const MKL_INT* ))

_Mkl_Api(int,vslsCorrExecX1D,(VSLCorrTaskPtr  , const float [],  const MKL_INT  , float [],  const MKL_INT  ))
_mkl_api(int,vslscorrexecx1d,(VSLCorrTaskPtr* , const float [],  const MKL_INT* , float [],  const MKL_INT* ))
_MKL_API(int,VSLSCORREXECX1D,(VSLCorrTaskPtr* , const float [],  const MKL_INT* , float [],  const MKL_INT* ))

_Mkl_Api(int,vslzCorrExecX1D,(VSLCorrTaskPtr  , const MKL_Complex16 [], const MKL_INT  , MKL_Complex16 [], const MKL_INT  ))
_mkl_api(int,vslzcorrexecx1d,(VSLCorrTaskPtr* , const MKL_Complex16 [], const MKL_INT* , MKL_Complex16 [], const MKL_INT* ))
_MKL_API(int,VSLZCORREXECX1D,(VSLCorrTaskPtr* , const MKL_Complex16 [], const MKL_INT* , MKL_Complex16 [], const MKL_INT* ))

_Mkl_Api(int,vslcCorrExecX1D,(VSLCorrTaskPtr  , const MKL_Complex8 [],  const MKL_INT  , MKL_Complex8 [],  const MKL_INT  ))
_mkl_api(int,vslccorrexecx1d,(VSLCorrTaskPtr* , const MKL_Complex8 [],  const MKL_INT* , MKL_Complex8 [],  const MKL_INT* ))
_MKL_API(int,VSLCCORREXECX1D,(VSLCorrTaskPtr* , const MKL_Complex8 [],  const MKL_INT* , MKL_Complex8 [],  const MKL_INT* ))


/*
//++
//  SUMMARARY STATTISTICS LIBRARY ROUTINES
//--
*/

/*
//  Task constructors
*/
_Mkl_Api(int,vsldSSNewTask,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const double [], const double [], const MKL_INT []))
_mkl_api(int,vsldssnewtask,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const double [], const double [], const MKL_INT []))
_MKL_API(int,VSLDSSNEWTASK,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const double [], const double [], const MKL_INT []))

_Mkl_Api(int,vslsSSNewTask,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const float  [], const float  [], const MKL_INT []))
_mkl_api(int,vslsssnewtask,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const float  [], const float  [], const MKL_INT []))
_MKL_API(int,VSLSSSNEWTASK,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* , const MKL_INT* , const float  [], const float  [], const MKL_INT []))


/*
// Task editors
*/

/*
// Editor to modify a task parameter
*/
_Mkl_Api(int,vsldSSEditTask,(VSLSSTaskPtr  , const MKL_INT  , const double* ))
_mkl_api(int,vsldssedittask,(VSLSSTaskPtr* , const MKL_INT* , const double* ))
_MKL_API(int,VSLDSSEDITTASK,(VSLSSTaskPtr* , const MKL_INT* , const double* ))

_Mkl_Api(int,vslsSSEditTask,(VSLSSTaskPtr  , const MKL_INT  , const float* ))
_mkl_api(int,vslsssedittask,(VSLSSTaskPtr* , const MKL_INT* , const float* ))
_MKL_API(int,VSLSSSEDITTASK,(VSLSSTaskPtr* , const MKL_INT* , const float* ))

_Mkl_Api(int,vsliSSEditTask,(VSLSSTaskPtr  , const MKL_INT  , const MKL_INT* ))
_mkl_api(int,vslissedittask,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* ))
_MKL_API(int,VSLISSEDITTASK,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* ))

/*
// Task specific editors
*/

/*
// Editors to modify moments related parameters
*/
_Mkl_Api(int,vsldSSEditMoments,(VSLSSTaskPtr  , const double* , const double* , const double* , const double* , const double* , const double* , const double* ))
_mkl_api(int,vsldsseditmoments,(VSLSSTaskPtr* , const double* , const double* , const double* , const double* , const double* , const double* , const double* ))
_MKL_API(int,VSLDSSEDITMOMENTS,(VSLSSTaskPtr* , const double* , const double* , const double* , const double* , const double* , const double* , const double* ))

_Mkl_Api(int,vslsSSEditMoments,(VSLSSTaskPtr  , const float* , const float* , const float* , const float* , const float* , const float* , const float* ))
_mkl_api(int,vslssseditmoments,(VSLSSTaskPtr* , const float* , const float* , const float* , const float* , const float* , const float* , const float* ))
_MKL_API(int,VSLSSSEDITMOMENTS,(VSLSSTaskPtr* , const float* , const float* , const float* , const float* , const float* , const float* , const float* ))


/*
// Editors to modify sums related parameters
*/
_Mkl_Api(int,vsldSSEditSums,(VSLSSTaskPtr  , const double* , const double* , const double* , const double* , const double* , const double* , const double* ))
_mkl_api(int,vsldsseditsums,(VSLSSTaskPtr* , const double* , const double* , const double* , const double* , const double* , const double* , const double* ))
_MKL_API(int,VSLDSSEDITSUMS,(VSLSSTaskPtr* , const double* , const double* , const double* , const double* , const double* , const double* , const double* ))

_Mkl_Api(int,vslsSSEditSums,(VSLSSTaskPtr  , const float* , const float* , const float* , const float* , const float* , const float* , const float* ))
_mkl_api(int,vslssseditsums,(VSLSSTaskPtr* , const float* , const float* , const float* , const float* , const float* , const float* , const float* ))
_MKL_API(int,VSLSSSEDITSUMS,(VSLSSTaskPtr* , const float* , const float* , const float* , const float* , const float* , const float* , const float* ))


/*
// Editors to modify variance-covariance/correlation matrix related parameters
*/
_Mkl_Api(int,vsldSSEditCovCor,(VSLSSTaskPtr  , const double* , const double* ,  const MKL_INT* , const double* , const MKL_INT* ))
_mkl_api(int,vsldsseditcovcor,(VSLSSTaskPtr* , const double* , const double* ,  const MKL_INT* , const double* , const MKL_INT* ))
_MKL_API(int,VSLDSSEDITCOVCOR,(VSLSSTaskPtr* , const double* , const double* ,  const MKL_INT* , const double* , const MKL_INT* ))

_Mkl_Api(int,vslsSSEditCovCor,(VSLSSTaskPtr  , const float* , const float* , const MKL_INT* , const float* , const MKL_INT* ))
_mkl_api(int,vslssseditcovcor,(VSLSSTaskPtr* , const float* , const float* , const MKL_INT* , const float* , const MKL_INT* ))
_MKL_API(int,VSLSSSEDITCOVCOR,(VSLSSTaskPtr* , const float* , const float* , const MKL_INT* , const float* , const MKL_INT* ))


/*
// Editors to modify cross-product matrix related parameters
*/
_Mkl_Api(int,vsldSSEditCP,(VSLSSTaskPtr  , const double* , const double* ,  const double* , const MKL_INT* ))
_mkl_api(int,vsldsseditcp,(VSLSSTaskPtr* , const double* , const double* ,  const double* , const MKL_INT* ))
_MKL_API(int,VSLDSSEDITCP,(VSLSSTaskPtr* , const double* , const double* ,  const double* , const MKL_INT* ))

_Mkl_Api(int,vslsSSEditCP,(VSLSSTaskPtr  , const float* , const float* , const float* , const MKL_INT* ))
_mkl_api(int,vslssseditcp,(VSLSSTaskPtr* , const float* , const float* , const float* , const MKL_INT* ))
_MKL_API(int,VSLSSSEDITCP,(VSLSSTaskPtr* , const float* , const float* , const float* , const MKL_INT* ))


/*
// Editors to modify partial variance-covariance matrix related parameters
*/
_Mkl_Api(int,vsldSSEditPartialCovCor,(VSLSSTaskPtr  , const MKL_INT [], const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* ))
_mkl_api(int,vsldsseditpartialcovcor,(VSLSSTaskPtr* , const MKL_INT [], const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* ))
_MKL_API(int,VSLDSSEDITPARTIALCOVCOR,(VSLSSTaskPtr* , const MKL_INT [], const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* ))

_Mkl_Api(int,vslsSSEditPartialCovCor,(VSLSSTaskPtr  , const MKL_INT [], const float* , const MKL_INT* , const float* , const MKL_INT* , const float* ,  const MKL_INT* , const float* ,  const MKL_INT* ))
_mkl_api(int,vslssseditpartialcovcor,(VSLSSTaskPtr* , const MKL_INT [], const float* , const MKL_INT* , const float* , const MKL_INT* , const float* ,  const MKL_INT* , const float* ,  const MKL_INT* ))
_MKL_API(int,VSLSSSEDITPARTIALCOVCOR,(VSLSSTaskPtr* , const MKL_INT [], const float* , const MKL_INT* , const float* , const MKL_INT* , const float* ,  const MKL_INT* , const float* ,  const MKL_INT* ))


/*
// Editors to modify quantiles related parameters
*/
_Mkl_Api(int,vsldSSEditQuantiles,(VSLSSTaskPtr  , const MKL_INT* , const double* , const double* , const double* , const MKL_INT* ))
_mkl_api(int,vsldsseditquantiles,(VSLSSTaskPtr* , const MKL_INT* , const double* , const double* , const double* , const MKL_INT* ))
_MKL_API(int,VSLDSSEDITQUANTILES,(VSLSSTaskPtr* , const MKL_INT* , const double* , const double* , const double* , const MKL_INT* ))

_Mkl_Api(int,vslsSSEditQuantiles,(VSLSSTaskPtr  , const MKL_INT* , const float* , const float* , const float* , const MKL_INT* ))
_mkl_api(int,vslssseditquantiles,(VSLSSTaskPtr* , const MKL_INT* , const float* , const float* , const float* , const MKL_INT* ))
_MKL_API(int,VSLSSSEDITQUANTILES,(VSLSSTaskPtr* , const MKL_INT* , const float* , const float* , const float* , const MKL_INT* ))


/*
// Editors to modify stream data quantiles related parameters
*/
_Mkl_Api(int,vsldSSEditStreamQuantiles,(VSLSSTaskPtr  , const MKL_INT* , const double* , const double* , const MKL_INT* , const double* ))
_mkl_api(int,vsldsseditstreamquantiles,(VSLSSTaskPtr* , const MKL_INT* , const double* , const double* , const MKL_INT* , const double* ))
_MKL_API(int,VSLDSSEDITSTREAMQUANTILES,(VSLSSTaskPtr* , const MKL_INT* , const double* , const double* , const MKL_INT* , const double* ))

_Mkl_Api(int,vslsSSEditStreamQuantiles,(VSLSSTaskPtr  , const MKL_INT* , const float* , const float* , const MKL_INT* , const float* ))
_mkl_api(int,vslssseditstreamquantiles,(VSLSSTaskPtr* , const MKL_INT* , const float* , const float* , const MKL_INT* , const float* ))
_MKL_API(int,VSLSSSEDITSTREAMQUANTILES,(VSLSSTaskPtr* , const MKL_INT* , const float* , const float* , const MKL_INT* , const float* ))

/*
// Editors to modify pooled/group variance-covariance matrix related parameters
*/
_Mkl_Api(int,vsldSSEditPooledCovariance,(VSLSSTaskPtr  , const MKL_INT* , const double* , const double* , const MKL_INT* , const double* , const double* ))
_mkl_api(int,vsldsseditpooledcovariance,(VSLSSTaskPtr* , const MKL_INT* , const double* , const double* , const MKL_INT* , const double* , const double* ))
_MKL_API(int,VSLDSSEDITPOOLEDCOVARIANCE,(VSLSSTaskPtr* , const MKL_INT* , const double* , const double* , const MKL_INT* , const double* , const double* ))

_Mkl_Api(int,vslsSSEditPooledCovariance,(VSLSSTaskPtr  , const MKL_INT* , const float* , const float* , const MKL_INT* , const float* , const float* ))
_mkl_api(int,vslssseditpooledcovariance,(VSLSSTaskPtr* , const MKL_INT* , const float* , const float* , const MKL_INT* , const float* , const float* ))
_MKL_API(int,VSLSSSEDITPOOLEDCOVARIANCE,(VSLSSTaskPtr* , const MKL_INT* , const float* , const float* , const MKL_INT* , const float* , const float* ))


/*
// Editors to modify robust variance-covariance matrix related parameters
*/
_Mkl_Api(int,vsldSSEditRobustCovariance,(VSLSSTaskPtr  , const MKL_INT* , const MKL_INT* ,  const double* , const double* , const double* ))
_mkl_api(int,vsldsseditrobustcovariance,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* ,  const double* , const double* , const double* ))
_MKL_API(int,VSLDSSEDITROBUSTCOVARIANCE,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* ,  const double* , const double* , const double* ))

_Mkl_Api(int,vslsSSEditRobustCovariance,(VSLSSTaskPtr  , const MKL_INT* , const MKL_INT* ,  const float* , const float* , const float* ))
_mkl_api(int,vslssseditrobustcovariance,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* ,  const float* , const float* , const float* ))
_MKL_API(int,VSLSSSEDITROBUSTCOVARIANCE,(VSLSSTaskPtr* , const MKL_INT* , const MKL_INT* ,  const float* , const float* , const float* ))


/*
// Editors to modify outliers detection parameters
*/
_Mkl_Api(int,vsldSSEditOutliersDetection,(VSLSSTaskPtr  , const MKL_INT* , const double* , const double* ))
_mkl_api(int,vsldsseditoutliersdetection,(VSLSSTaskPtr* , const MKL_INT* , const double* , const double* ))
_MKL_API(int,VSLDSSEDITOUTLIERSDETECTION,(VSLSSTaskPtr* , const MKL_INT* , const double* , const double* ))

_Mkl_Api(int,vslsSSEditOutliersDetection,(VSLSSTaskPtr  , const MKL_INT* , const float* , const float* ))
_mkl_api(int,vslssseditoutliersdetection,(VSLSSTaskPtr* , const MKL_INT* , const float* , const float* ))
_MKL_API(int,VSLSSSEDITOUTLIERSDETECTION,(VSLSSTaskPtr* , const MKL_INT* , const float* , const float* ))

/*
// Editors to modify missing values support parameters
*/
_Mkl_Api(int,vsldSSEditMissingValues,(VSLSSTaskPtr  , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* ))
_mkl_api(int,vsldsseditmissingvalues,(VSLSSTaskPtr* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* ))
_MKL_API(int,VSLDSSEDITMISSINGVALUES,(VSLSSTaskPtr* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* , const MKL_INT* , const double* ))

_Mkl_Api(int,vslsSSEditMissingValues,(VSLSSTaskPtr  , const MKL_INT* , const float* , const MKL_INT* , const float* , const MKL_INT* , const float* , const MKL_INT* , const float* , const MKL_INT* , const float* ))
_mkl_api(int,vslssseditmissingvalues,(VSLSSTaskPtr* , const MKL_INT* , const float* , const MKL_INT* , const float* , const MKL_INT* , const float* , const MKL_INT* , const float* , const MKL_INT* , const float* ))
_MKL_API(int,VSLSSSEDITMISSINGVALUES,(VSLSSTaskPtr* , const MKL_INT* , const float* , const MKL_INT* , const float* , const MKL_INT* , const float* , const MKL_INT* , const float* , const MKL_INT* , const float* ))

/*
// Editors to modify matrixparametrization parameters
*/
_Mkl_Api(int,vsldSSEditCorParameterization,(VSLSSTaskPtr  , const double* , const MKL_INT* , const double* , const MKL_INT* ))
_mkl_api(int,vsldsseditcorparameterization,(VSLSSTaskPtr* , const double* , const MKL_INT* , const double* , const MKL_INT* ))
_MKL_API(int,VSLDSSEDITCORPARAMETERIZATION,(VSLSSTaskPtr* , const double* , const MKL_INT* , const double* , const MKL_INT* ))

_Mkl_Api(int,vslsSSEditCorParameterization,(VSLSSTaskPtr  , const float* , const MKL_INT* , const float* , const MKL_INT* ))
_mkl_api(int,vslssseditcorparameterization,(VSLSSTaskPtr* , const float* , const MKL_INT* , const float* , const MKL_INT* ))
_MKL_API(int,VSLSSSEDITCORPARAMETERIZATION,(VSLSSTaskPtr* , const float* , const MKL_INT* , const float* , const MKL_INT* ))


/*
// Compute routines
*/
_Mkl_Api(int,vsldSSCompute,(VSLSSTaskPtr  , const unsigned MKL_INT64  , const MKL_INT  ))
_mkl_api(int,vsldsscompute,(VSLSSTaskPtr* , const unsigned MKL_INT64* , const MKL_INT* ))
_MKL_API(int,VSLDSSCOMPUTE,(VSLSSTaskPtr* , const unsigned MKL_INT64* , const MKL_INT* ))

_Mkl_Api(int,vslsSSCompute,(VSLSSTaskPtr  , const unsigned MKL_INT64  , const MKL_INT  ))
_mkl_api(int,vslssscompute,(VSLSSTaskPtr* , const unsigned MKL_INT64* , const MKL_INT* ))
_MKL_API(int,VSLSSSCOMPUTE,(VSLSSTaskPtr* , const unsigned MKL_INT64* , const MKL_INT* ))


/*
// Task destructor
*/
_Mkl_Api(int,vslSSDeleteTask,(VSLSSTaskPtr* ))
_mkl_api(int,vslssdeletetask,(VSLSSTaskPtr* ))
_MKL_API(int,VSLSSDELETETASK,(VSLSSTaskPtr* ))

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __MKL_VSL_FUNCTIONS_H__ */
