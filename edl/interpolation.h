#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "edl/edl.h"
#include "edl/edlerror.h"

namespace EDL_NAMESPACE
{

template <typename VALUE_TYPE>
class LinearInterpolationFunction
{

private:

  real       m_X[2];
  VALUE_TYPE m_Y[2];


public:

  VALUE_TYPE compute(real x)
  {
    real w = (x - m_X[0])/(m_X[1] - m_X[0]);
    return (1-w)*m_Y[0] + w*m_Y[1];
  }

  int size() { return 2; }
  int offset() { return 0; }
  void update() {}

  void set(size_t i, real x, VALUE_TYPE y)
  {
    m_X[i] = x;
    m_Y[i] = y;
  }

};

template <typename STENCIL_FUNCTION, typename VALUE_TYPE>
class Interpolation
{

private: // attributes

  TList<real>       *m_XList;
  TList<VALUE_TYPE> *m_YList;
  STENCIL_FUNCTION   m_Stencil;
  real               m_LowerBound;
  real               m_UpperBound;
  size_t             m_Idx;
  bool               m_FirstCall;
  bool               m_AllowExtrapolation;


private: // methods

  void findIndex(real x)
  {
    bool search_down = false;
    if (m_FirstCall) {
      m_FirstCall = false;
      m_Idx = 0;
    } else {
      if (x < m_Stencil.m_X[0]) {
        search_down = true;
        m_Idx = m_XList->prevIdx(m_Idx);
      } else {
        m_Idx = m_XList->nextIdx(m_Idx);
      }
    }
    bool found = false;
    if (search_down) {
      while (m_Idx >= m_XList->beginIdx()) {
        if (m_XList->at(m_Idx) <= x) {
          found = true;
          break;
        }
        m_Idx = m_XList->prevIdx(m_Idx);
      }
    } else {
      while (m_Idx < m_XList->lastIdx()) {
        size_t idx = m_XList->nextIdx(m_Idx);
        if (m_XList->at(idx) >= x) {
          found = true;
          break;
        }
        m_Idx = idx;
      }
    }
    if (!found && !m_AllowExtrapolation) {
      throw EdlError("x value out of bounds");
    }
    size_t idx = m_Idx;
    for (int i = 0; i < m_Stencil.offset(); ++i) {
      if (idx > m_XList->beginIdx()) {
        idx = m_XList->prevIdx(idx);
      }
    }
    for (int i = 0; i < m_Stencil.size(); ++i) {
      if (idx < m_XList->endIdx()) {
        m_Stencil.set(i, m_XList->at(idx), m_YList->qt(idx));
        idx = m_XList->nexIdx(idx);
      } else {
        throw EdlError("index out of bounds");
      }
    }
  }


public:

  Interpolator(TList<real> *x_list, TList<VALUE_TYPE> *y_list, bool allow_extrapolation = false)
  {
    m_XList = x_list;
    m_YList = y_list;
    m_FirstCall = true;
    m_AllowExtrapolation = allow_extrapolation;
  }

  VALUE_TYPE operator()(real x)
  {
    if (m_FirstCall) {
      findIndex(x);
    } else if (x >= m_LowerBound && x <= m_UpperBound) {
      return compute(x);
    } else {
      findIndex(x);
    }
    return m_Stencil.compute(x);
  }

};

}

#endif // INTERPOLATION_H
