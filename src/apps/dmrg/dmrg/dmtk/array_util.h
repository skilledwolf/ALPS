/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2006 -2010 by Adrian Feiguin <afeiguin@uwyo.edu>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef __DMTK_ARRAY_UTIL_H__
#define __DMTK_ARRAY_UTIL_H__

#include <complex>
#include "conj.h"
#ifdef WITH_LAPACK
#include "lapack_interface.h"
#endif

#ifdef WITH_CBLAS
extern "C" {
#include <atlas_misc.h>
#include <atlas_cblastypealias.h>
#include <atlas_level1.h>
#include <atlas_level2.h>
#include <atlas_level3.h>
#include <cblas.h>
#ifdef WITH_PTHREADS
   #include "atlas_ptalias1.h"
   #include "atlas_ptalias2.h"
   #include "atlas_ptalias3.h"
#endif
}
#endif // WITH_CBLAS

using namespace std;

namespace dmtk
{

#ifdef WITH_LAPACK

double
dot_product(const int& n, 
            const double *v1, const int& stride1, 
            const double *v2, const int& stride2)
{
#ifdef WITH_CBLAS
  return cblas_ddot(n, v1, stride1, v2, stride2);  
#else
  return FORTRAN_ID(ddot)(n, v1, stride1, v2, stride2);  
#endif
}

complex<double>
dot_product(const int& n, 
            const complex<double> *v1, const int& stride1, 
            const complex<double> *v2, const int& stride2)
{
  complex<double> sum(0.f,0.f);
  const complex<double> *pi = v1;
  const complex<double> *pj = v2;

//  sum = FORTRAN_ID(zdotc)(n, v1, stride1, v2, stride2);  
//  return sum;

  if(stride1 == 1 && stride2 == 1)
    for(int k = 0; k < n; k++)
      sum += std::conj(*pi++)*(*pj++);
  else
    for(int k = 0; k < n; k++){
      sum += std::conj(*pi)*(*pj);
      pi += stride1;
      pj += stride2;
    }

  return sum;
}

void
array_copy(int n, const double *in, double *out)
{
#ifdef WITH_CBLAS
  cblas_dcopy(n, in, 1, out, 1);
#else
  FORTRAN_ID(dcopy)(n, in, 1, out, 1);
#endif
}

void
array_copy(int n, const complex<double> *in, complex<double> *out)
{
#ifdef WITH_CBLAS
  cblas_zcopy(n, in, 1, out, 1);
#else
  FORTRAN_ID(zcopy)(n, in, 1, out, 1);
#endif
}


#else // !WITH_LAPACK

template<class T>
T
dot_product(const int& n, 
            const T *v1, const int& stride1, const T *v2, const int& stride2)
{
  T sum(0);
  const T *pi = v1;
  const T *pj = v2;

  if(stride1 == 1 && stride2 == 1)
    for(int k = 0; k < n; k++)
      sum += std::conj(*pi++)*(*pj++);
  else
    for(int k = 0; k < n; k++){
      sum += std::conj(*pi)*(*pj);
      pi += stride1;
      pj += stride2;
    }

  return sum;
}


#endif // WITH_LAPACK

template<class T>
void array_copy(int n, const T* in, T* out)
{
  for (int i = 0; i < n; ++i) out[i] = in[i];
}


template<class T>
void
array_copy(int n, T& in, T& out)
{
  for (int i = 0; i < n; ++i) out[i] = in[i];
}

template<class T1, class T2>
void
array_copy2(int n, const T1& in, T2& out)
{
  for (int i = 0; i < n; ++i) out[i] = in[i];
}

template<class T1, class T2>
void
array_copy2(int n, T1& in, T2& out)
{
  for (int i = 0; i < n; ++i) out[i] = in[i];
}

inline double quickran(long & idum)
{
    const int im = 134456;
    const int ia = 8121;
    const int ic = 28411;
    const double scale = 1.0 / im;
    idum = (idum*ia+ic)%im;
    return double(idum) * scale;
}

} // namespace dmtk

#endif // __DMTK_ARRAY_UTIL_H__
