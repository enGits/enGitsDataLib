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
#include "interpolationcurve.h"

namespace EDL_NAMESPACE
{

InterpolationCurve::InterpolationCurve()
{
  m_SubSamplingSteps = 0;
  m_InterpolationOrder = 3;
  m_OriginalY.link(&m_OriginalX);
  m_InterY.link(&m_InterX);
  m_Interp = NULL;
  addPoint(0, 1);
  addPoint(1, 1);
  update();
}

InterpolationCurve::InterpolationCurve(const InterpolationCurve &other_curve)
{
  if (!other_curve.m_OriginalX.isClean()) {
    throw EdlError("trying to copy a non-clean InterpolationCurve");
  }
  m_SubSamplingSteps   = other_curve.m_SubSamplingSteps;
  m_InterpolationOrder = other_curve.m_InterpolationOrder;
  m_OriginalY.link(&m_OriginalX);
  m_InterY.link(&m_InterX);
  for (size_t i = 0; i < other_curve.m_OriginalX.numEntries(); ++i) {
    size_t i_new = m_OriginalX.addEntry();
    m_OriginalX[i_new] = other_curve.m_OriginalX[i];
    m_OriginalY[i_new] = other_curve.m_OriginalY[i];
  }
  m_Interp = NULL;
  update();
}

void InterpolationCurve::addPoint(real x, real y)
{
  size_t i = m_OriginalX.addEntry();
  m_OriginalX[i] = x;
  m_OriginalY[i] = y;
}

void InterpolationCurve::update()
{
  if (m_OriginalX.numActiveEntries() < 2) {
    return;
  }

  // compute derivatives
  //
  TList<real> grad(&m_OriginalX);
  m_OriginalX.cleanUp();
  for (size_t i = 0; i < m_OriginalX.numEntries(); ++i) {
    size_t i1 = i;
    size_t i2 = i + 1;
    if (i > 0) {
      --i1;
    }
    if (i == m_OriginalX.numEntries() - 1) {
      --i2;
    }
    grad[i] = (m_OriginalY[i2] - m_OriginalY[i1])/(m_OriginalX[i2] - m_OriginalX[i1]);
  }
  grad[0] = 2*grad[0] - grad[1];
  grad[m_OriginalX.numEntries() - 1] = 2*grad[m_OriginalX.numEntries() - 1] - grad[m_OriginalX.numEntries() - 2];

  // compute sub-sampling with Bezier curves
  //
  m_InterX.delAll();
  if (m_InterpolationOrder < 1 || m_InterpolationOrder > 3) {
    throw EdlError("unsupported interpolation order");
  }
  if (m_InterpolationOrder > 1) {
    for (size_t i = 0; i < m_OriginalX.numEntries() - 1; ++i) {
      size_t i1 = i;
      size_t i2 = i + 1;
      real x1 = m_OriginalX[i1];
      real y1 = m_OriginalY[i1];
      real x4 = m_OriginalX[i2];
      real y4 = m_OriginalY[i2];
      real x2 = 0.5*(x1 + x4);
      real y2 = 0.5*(y1 + y4);
      real x3 = x2;
      real y3 = y2;
      real m1 = grad[i1];
      real m4 = grad[i2];
      if (fabs(m1 - m4) > 1e-6*(fabs(m1) + fabs(m4))) {
        if (m_InterpolationOrder == 2) {
          x2 = (y4 + m1*x1 - y1 - m4*x4)/(m1 - m4);
          y2 = 0.5*(y1 + m1*(x2 - x1) + y4 + m4*(x2 - x4));
        } else {
          real w = 0.25;
          x2 = (1-w)*x1 +    w *x4;
          x3 =    w *x1 + (1-w)*x4;
          y2 = y1 + w*(x4-x1)*m1;
          y3 = y4 - w*(x4-x1)*m4;
        }
      }
      int N = std::max(1, m_SubSamplingSteps + 1);
      real dw = 1.0/N;
      for (int j = 0; j < N; ++j) {
        size_t k = m_InterX.addEntry();
        real w = j*dw;
        real x = cube(1-w)*x1 + 3*sqr(1-w)*w*x2 + 3*(1-w)*sqr(w)*x3 + cube(w)*x4;
        real y = cube(1-w)*y1 + 3*sqr(1-w)*w*y2 + 3*(1-w)*sqr(w)*y3 + cube(w)*y4;
        if (x < x1) {
          throw EdlError("error while constructing interpolated data");
        }
        m_InterX[k] = x;
        m_InterY[k] = y;
      }
    }
    {
      size_t k = m_InterX.addEntry();
      m_InterX[k] = m_OriginalX.lastEntry();
      m_InterY[k] = m_OriginalY.lastEntry();
    }
  } else {
    FORALL(i, m_OriginalX.) {
      size_t j = m_InterX.addEntry();
      m_InterX[j] = m_OriginalX[i];
      m_InterY[j] = m_OriginalY[i];
    }
  }

  delete m_Interp;
  m_Interp = new Interpolation<LinInterpFunc<real> >(&m_InterX, &m_InterY, true);

  // m_OriginalX.delAll();
}

real InterpolationCurve::interpolate(real x)
{
  return m_Interp->interpolate(x);
}

void InterpolationCurve::clear()
{
   m_OriginalX.delAll();
}

}

