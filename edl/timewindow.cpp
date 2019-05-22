// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2015-2019 enGits GmbH                                    +
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

#include "timewindow.h"

namespace EDL_NAMESPACE
{

TimeWindow::TimeWindow(edl::real dt, size_t max_num_entries, size_t delta_entries)
  : TQueue<edl::real>(max_num_entries, delta_entries)
{
  m_Dt       = dt;
  m_Integral = 0;
  //
  m_ValueList.link(this);
}

void TimeWindow::append(edl::real t0, edl::real v0)
{
  EDL_BUG;
  //
  //    X----------X----------X-- ... --X----------X----------X
  //    ^          ^                               ^          ^
  //    |          |                               |          |
  //    0          1                               2          3
  //
  int       i2 = TQueue<edl::real>::at(m_Oldest);
  edl::real t2 = TQueue<edl::real>::m_ValueList->at(i2);
  edl::real v2 = m_ValueList[i2];
  //
  int       i1 = m_Newest;
  edl::real t1 = TQueue<edl::real>::m_ValueList->at(i1);
  edl::real v1 = m_ValueList[i1];
  int       i0 = TQueue<edl::real>::insert(t0);
  //
  m_ValueList[i0] = v0;
  //
  m_Integral += 0.5*(t0 - t1)*(v0 + v1);
  //
  if (t0 - t2 > m_Dt) {
    int       i3 = m_Oldest;
    edl::real t3 = TQueue<edl::real>::m_ValueList->at(i3);
    edl::real v3 = m_ValueList[i3];
    //
    m_Integral -= 0.5*(t2 - t3)*(v2 + v3);
  }
}

}
