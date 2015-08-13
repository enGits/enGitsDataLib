// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2015 enGits GmbH                                         +
// +                                                                    +
// + enGitsDataLib is free software: you can redistribute it and/or     +
// + modify it under the terms of the GNU Lesser General Public License +
// + as published by the Free Software Foundation, either version 3 of  +
// + the License, or (at your option) any later version.                +
// +                                                                    +
// + enGitsDataLib is distributed in the hope that it will be useful,   +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of     +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      +
// + GNU Lesser General Public License for more details.                +
// +                                                                    +
// + You should have received a copy of the GNU Lesser General Public   +
// + License along with enGitsDataLib.                                  +
// + If not, see <http://www.gnu.org/licenses/>.                        +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef EDL_H
#define EDL_H

#define EDL_NAMESPACE edl

#include <cfloat>

namespace EDL_NAMESPACE
{
#ifdef REAL_FLOAT
  typedef float real;
  const real max_real   =  FTL_MAX;
  const real min_real   = -FTL_MAX;
  const real small_real =  FLT_MIN;
#else
  typedef double real;
  const real max_real   =  DBL_MAX;
  const real min_real   = -DBL_MAX;
  const real small_real =  DBL_MIN;
#endif

  inline real sign(real x)
  {
    if (x >= 0) return  1;
    return -1;
  }

  inline real sign0(real x)
  {
    if      (x > 0) return  1;
    else if (x < 0) return -1;
    return 0;
  }

  inline real sqr(real x)
  {
    return x*x;
  }

  inline real ipow(real x, int e)
  {
    real v = 1;
    for (int i = 0; i < e; ++i) {
      v *= x;
    }
    return v;
  }

  inline real realmod(real a, real b)
  {
    int c = int(a/b);
    return a - c*b;
  }

}

#endif
