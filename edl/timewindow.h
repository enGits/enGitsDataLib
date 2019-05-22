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

#ifndef TIMEWINDOW_H
#define TIMEWINDOW_H

class TimeWindow;

#include "edl/tqueue.h"

namespace EDL_NAMESPACE
{

class TimeWindow : public TQueue<edl::real>
{

protected: // attributes

  edl::real             m_Dt;
  edl::real             m_Integral;
  edl::TList<edl::real> m_ValueList;

public: // methods

  TimeWindow(edl::real dt, size_t max_num_entries=100, size_t delta_entries=100);

  void append(edl::real t, edl::real value);

  edl::real average() { return m_Integral/m_Dt; }

};

}

#endif
