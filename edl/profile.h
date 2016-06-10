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

#include "edl.h"
#ifdef WITH_CXX11
#include <chrono>
#endif

namespace EDL_NAMESPACE
{

class Profile
{

#ifdef WITH_CXX11
private:

  typedef std::chrono::steady_clock clock_t;

  static clock_t::time_point m_Time[10];
  static clock_t::duration   m_Duration[10];


public:

  static int    size()       { return 10; }
  static double time(int i)  { return std::chrono::duration_cast<std::chrono::duration<double> >(m_Duration[i]).count(); }
  static void   start(int i) { m_Time[i] = clock_t::now(); }
  static void   stop(int i)  { m_Duration[i] += (clock_t::now() - m_Time[i]); }
#endif

};

}
