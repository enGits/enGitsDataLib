// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef INTERPOLATIONCURVE_H
#define INTERPOLATIONCURVE_H

#include "edl/interpolation.h"
#include "edl/tlist.h"

namespace EDL_NAMESPACE
{

class InterpolationCurve
{

private: // attributes

  TTList<real,100,100>                 m_OriginalX;
  TTList<real,100,100>                 m_OriginalY;
  TTList<real,1001,1000>               m_InterX;
  TTList<real,1001,1000>               m_InterY;
  Interpolation<LinInterpFunc<real> > *m_Interp;
  int                                  m_SubSamplingSteps;
  int                                  m_InterpolationOrder;


private: // methods

  edl::real sqr(edl::real x) { return x*x; }
  edl::real cube(edl::real x) { return x*x*x; }


public: // methods

  InterpolationCurve();
  InterpolationCurve(const InterpolationCurve& other_curve);
  void addPoint(real x, real y);
  void update();
  real interpolate(real x);
  real minX() { return m_InterX.firstEntry(); }
  real maxX() { return m_InterX.lastEntry(); }
  void clear();
  size_t numOriginalPoints() { return m_OriginalX.numEntries(); }
  size_t numInterpolatedPoints() { return m_InterX.numEntries(); }
  real getInterpolatedX(size_t i) { return m_InterX[i]; }
  real getInterpolatedY(size_t i) { return m_InterY[i]; }
  real getOriginalX(size_t i) { return m_OriginalX[i]; }
  real getOriginalY(size_t i) { return m_OriginalY[i]; }

};

}

#endif // INTERPOLATIONCURVE_H
