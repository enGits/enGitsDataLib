// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef EDL_H
#define EDL_H

#define EDL_NAMESPACE edl

namespace EDL_NAMESPACE
{
template <class V> struct MathVector;
template <class T, unsigned int DIM> class StaticVector;
}

#ifdef _WIN32
#undef min
#undef max
#endif

#include <cstdint>
#include <cfloat>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>

#include "edl/doctest.h"

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

  inline real snap(real v, real v0, real Dv)
  {
    if (fabs(v - v0) <= Dv) {
      return v0;
    }
    return v;
  }

  inline real lowLimit(real v, real v0)
  {
    if (v < v0) {
      return v0;
    }
    return v;
  }

  inline real highLimit(real v, real v0)
  {
    if (v > v0) {
      return v0;
    }
    return v;
  }

  inline real bandLimit(real v, real v1, real v2)
  {
    if (v < v1) {
      return v1;
    }
    if (v > v2) {
      return v2;
    }
    return v;
  }

  inline real linearSwitch(real x, real x1, real x2, real y1=-1.0, real y2=1.0)
  {
    real a = (y2 - y1)/(x2 - x1);
    real b = y1 - a*x1;
    real y_min, y_max;
    if (y1 > y2) {
      y_max = y1;
      y_min = y2;
    } else {
      y_min = y1;
      y_max = y2;
    }
    return std::max(y_min, std::min(y_max, a*x + b));
  }

  inline real when(bool condition, real true_value, real false_value)
  {
    if (condition) {
      return true_value;
    }
    return false_value;
  }

  struct smartBreakPoint
  {
    static std::map<std::string,uint64_t> all_hits;
    std::string name;
    uint64_t hits;
    smartBreakPoint(std::string name, int64_t min_count=std::numeric_limits<int64_t>::max()) : name(name)
    {
      if (all_hits.find(name) == all_hits.end()) {
        all_hits[name] = 0;
      }
      all_hits[name]++;
      hits = all_hits[name];
      if (all_hits[name] >= min_count) {
        conditionalBreak();
      }
    }
    void conditionalBreak()
    {
      // Put a debugger breakpoint here!
      [[maybe_unused]] int dummy = 0; // [[maybe_unused]] tells the compiler to not warn about unused variable
    }
  };

  template<typename T>
  bool almostEqual(T x, T y, int ulp=1000)
  {
    // Example code shamelessly stolen from:
    // https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
    //
    // The machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
      // unless the result is subnormal
      || std::abs(x-y) < std::numeric_limits<T>::min();
  }

  template<> inline
  bool almostEqual(bool x, bool y, int)
  {
    return x == y;
  }

  template<> inline
  bool almostEqual(int x, int y, int)
  {
    return x == y;
  }

  template<typename T, unsigned int N> inline
  bool almostEqual(MathVector<StaticVector<T,N>> x, MathVector<StaticVector<T,N>> y, int ulp=1000)
  {
    bool almost_equal = true;
    for (int i = 0; i < N; ++i) {
      if (!almostEqual(x[i], y[i], ulp)) {
        almost_equal = false;
        break;
      }
    }
    return almost_equal;
  }

  template <typename T>
  inline T absmin(T x, T x_min)
  {
    if (x > 0) return std::min(x,  x_min);
    else       return std::max(x, -x_min);
  }

  template <typename T>
  inline T absmax(T x, T x_max)
  {
    if (x > 0) return std::max(x,  x_max);
    else       return std::min(x, -x_max);
  }

  template <typename T>
  bool inInterval(T x, T x1, T x2, int ulp=1000)
  {
    if (x >= x1 && x <= x2) {
      return true;
    }
    if (almostEqual(x, x1)) {
      return true;
    }
    if (almostEqual(x, x2)) {
      return true;
    }
    return false;
  }

}

#if defined(_MSC_VER) && (_MSC_VER < 1700) // less than Visual Studio 2012
  inline int isnan(edl::real x)
  {
    return _isnan(x);
  };
  inline int isinf(edl::real x)
  {
    return !_finite(x);
  };
#endif

template <typename T>
inline std::uint8_t ilog2(T x)
{
  std::uint8_t n = 0;
  while (x >>= 1) ++n;
  return n;
}

#endif // EDL_H


// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

TEST_CASE("ilog2")
{
  using namespace EDL_NAMESPACE;
  //
  CHECK(ilog2(0) == 0);
  CHECK(ilog2(1) == 0);
  CHECK(ilog2(2) == 1);
  CHECK(ilog2(4) == 2);
  CHECK(ilog2(8) == 3);
  CHECK(ilog2(16) == 4);
  CHECK(ilog2(32) == 5);
  CHECK(ilog2(64) == 6);
  CHECK(ilog2(128) == 7);
  CHECK(ilog2(256) == 8);
  CHECK(ilog2(512) == 9);
  CHECK(ilog2(1024) == 10);
  CHECK(ilog2(2048) == 11);
  CHECK(ilog2(4096) == 12);
  CHECK(ilog2(8192) == 13);
  CHECK(ilog2(16384) == 14);
  CHECK(ilog2(32768) == 15);
  CHECK(ilog2(65536) == 16);
  //
  CHECK(ilog2(2+1) == 1);
  CHECK(ilog2(4+1) == 2);
  CHECK(ilog2(8+1) == 3);
  CHECK(ilog2(16+1) == 4);
  CHECK(ilog2(32+1) == 5);
  CHECK(ilog2(64+1) == 6);
  CHECK(ilog2(128+1) == 7);
  CHECK(ilog2(256+1) == 8);
  CHECK(ilog2(512+1) == 9);
  CHECK(ilog2(1024+1) == 10);
  CHECK(ilog2(2048+1) == 11);
  CHECK(ilog2(4096+1) == 12);
  CHECK(ilog2(8192+1) == 13);
  CHECK(ilog2(16384+1) == 14);
  CHECK(ilog2(32768+1) == 15);
  CHECK(ilog2(65536+1) == 16);

}