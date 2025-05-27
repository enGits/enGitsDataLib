// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
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
