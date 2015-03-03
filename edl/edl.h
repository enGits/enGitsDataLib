#ifndef EDL_H
#define EDL_H

#define EDL_NAMESPACE edl

#include <cfloat>

#define forceinline __inline__ __attribute__((always_inline))
#define ensure_forceinline __attribute__((always_inline)) // inline or die

namespace EDL_NAMESPACE
{
#ifdef REAL_FLOAT
  typedef float real;
  const real max_real =  FTL_MAX;
  const real min_real = -FTL_MAX;
#else
  typedef double real;
  const real max_real =  DBL_MAX;
  const real min_real = -DBL_MAX;
#endif
}

#endif
