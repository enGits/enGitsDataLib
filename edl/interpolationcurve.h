// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2015-2016 enGits GmbH                                    +
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
  void clear();

};

}

#endif // INTERPOLATIONCURVE_H
