// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "edl.h"
#if __cplusplus >= 201103L
#include <chrono>
#endif

namespace EDL_NAMESPACE
{

class Profile
{

#if __cplusplus >= 201103L
private:

  typedef std::chrono::steady_clock clock_t;

  static clock_t::time_point m_Time[10];
  static clock_t::duration   m_Duration[10];


public:

  static int    size()       { return 10; }
  static double time(int i)  { return std::chrono::duration_cast<std::chrono::duration<double> >(m_Duration[i]).count(); }
  static void   start(int i) { m_Time[i] = clock_t::now(); }
  static void   stop(int i)  { m_Duration[i] += (clock_t::now() - m_Time[i]); }
#else
public:

  static int    size()       { return -1; }
  static double time(int i)  { return -1; }
  static void   start(int i) { }
  static void   stop(int i)  { }
#endif

};

}
