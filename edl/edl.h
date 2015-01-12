#ifndef EDL_H
#define EDL_H

#define EDL_NAMESPACE edl

#include <cfloat>

namespace EDL_NAMESPACE
{
#ifdef REAL_FLOAT
  typedef float real;
  const real max_real = FTL_MAX;
  const real min_real = FTL_MIN;
#else
  typedef double real;
  const real max_real = DBL_MAX;
  const real min_real = DBL_MIN;
#endif
}

#endif
