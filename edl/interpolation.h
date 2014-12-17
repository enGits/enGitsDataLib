#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "edl/edl.h"
#include "edl/edlerror.h"

namespace EDL_NAMESPACE
{

template <typename VALUE_TYPE>
struct StencilFunction
{

  typedef VALUE_TYPE value_t;

  real    m_X[2];
  value_t m_Y[2];
  value_t m_DyDx[2];

  int size() { return 2; }

  edl::real sqr(edl::real x) { return x*x; }
  edl::real cube(edl::real x) { return x*x*x; }

  void set(size_t i, edl::real x, value_t y)
  {
    m_X[i] = x;
    m_Y[i] = y;
  }

  void update() {}

};


template <typename VALUE_TYPE>
struct Bezier0 : public StencilFunction<VALUE_TYPE>
{
  typedef typename StencilFunction<VALUE_TYPE>::value_t value_t;
  using StencilFunction<VALUE_TYPE>::m_X;
  using StencilFunction<VALUE_TYPE>::m_Y;
  using StencilFunction<VALUE_TYPE>::m_DyDx;

  value_t compute(real x)
  {
    real t = (x - m_X[0])/(m_X[1] - m_X[0]);
    if (t < 0.5) return m_Y[0];
    return m_Y[1];
  }

};

template <typename VALUE_TYPE>
struct Bezier1 : public StencilFunction<VALUE_TYPE>
{
  typedef typename StencilFunction<VALUE_TYPE>::value_t value_t;
  using StencilFunction<VALUE_TYPE>::m_X;
  using StencilFunction<VALUE_TYPE>::m_Y;
  using StencilFunction<VALUE_TYPE>::m_DyDx;

  value_t compute(real x)
  {
    real t = (x - m_X[0])/(m_X[1] - m_X[0]);
    return (1-t)*m_Y[0] + t*m_Y[1];
  }

};

template <typename VALUE_TYPE>
struct Bezier2 : public StencilFunction<VALUE_TYPE>
{
  typedef typename StencilFunction<VALUE_TYPE>::value_t value_t;
  using StencilFunction<VALUE_TYPE>::m_X;
  using StencilFunction<VALUE_TYPE>::m_Y;
  using StencilFunction<VALUE_TYPE>::m_DyDx;

  edl::real m_Ye1;

  value_t compute(real x)
  {
    real t = (x - m_X[0])/(m_X[1] - m_X[0]);
    return (1-t)*((1-t)*m_Y[0]+t*m_Ye1) + t*((1-t)*m_Ye1+t*m_Y[1]);
  }

  void update()
  {
    m_Ye1 = m_Y[0] + 0.5*(m_Y[0] + 0.5*(m_X[1] - m_X[0])*m_DyDx[0] + m_Y[1] - 0.5*(m_X[1] - m_X[0])*m_DyDx[1]);
  }

};


template <typename STENCIL_FUNCTION>
class Interpolation
{

public: // data types

  typedef typename STENCIL_FUNCTION::value_t value_t;



private: // attributes

  TList<real>       *m_XList;
  TList<value_t>    *m_YList;
  STENCIL_FUNCTION   m_Stencil;
  real               m_LowerBound;
  real               m_UpperBound;
  size_t             m_Idx;
  bool               m_FirstCall;
  bool               m_AllowExtrapolation;


private: // methods

  void findIndex(real x)
  {
    if (!m_XList->isClean()) {
      m_XList->cleanUp();
    }
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
      while (m_Idx >= 0) {
        if (m_XList->at(m_Idx) <= x) {
          found = true;
          break;
        }
        --m_Idx;
      }
    } else {
      while (m_Idx < m_XList->lastIdx()) {
        size_t idx = m_Idx + 1;
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
    for (int i = 0; i < 2; ++i) {
      if (idx < m_XList->endIdx()) {
        m_Stencil.set(i, m_XList->at(idx), m_YList->at(idx));
        idx = m_XList->nextIdx(idx);
      } else {
        throw EdlError("index out of bounds");
      }

      if (m_Idx < 0 || m_Idx >= m_XList->lastIdx()) {
        throw EdlError("index out of bounds");
      }
      m_Stencil.m_X[0] = m_XList->at(m_Idx);
      m_Stencil.m_X[1] = m_XList->at(m_Idx + 1);
      m_Stencil.m_Y[0] = m_YList->at(m_Idx);
      m_Stencil.m_Y[1] = m_YList->at(m_Idx + 1);
      m_Stencil.m_DyDx[0] = 1.0/(m_Stencil.m_X[1] - m_Stencil.m_X[0])*(m_Stencil.m_Y[1] - m_Stencil.m_Y[0]);
      m_Stencil.m_DyDx[1] = m_Stencil.m_DyDx[0];
      if (m_Idx > 0) {
        m_Stencil.m_DyDx[0] = 1.0/(m_XList->at(m_Idx + 1) - m_XList->at(m_Idx - 1))*(m_YList->at(m_Idx + 1) - m_YList->at(m_Idx - 1));
      }
      if (m_Idx < m_XList->lastIdx() - 1) {
        m_Stencil.m_DyDx[1] = 1.0/(m_XList->at(m_Idx + 2) - m_XList->at(m_Idx))*(m_YList->at(m_Idx + 2) - m_YList->at(m_Idx));
      }

    }
    m_Stencil.update();
    m_LowerBound = m_Stencil.m_X[0];
    m_UpperBound = m_Stencil.m_X[1];
  }


public:

  Interpolation(TList<real> *x_list, TList<value_t> *y_list, bool allow_extrapolation = false)
  {
    m_XList = x_list;
    m_YList = y_list;
    m_FirstCall = true;
    m_AllowExtrapolation = allow_extrapolation;
  }

  value_t interpolate(real x)
  {
    if (m_FirstCall) {
      findIndex(x);
    } else if (x >= m_LowerBound && x <= m_UpperBound) {
      return m_Stencil.compute(x);
    } else {
      findIndex(x);
    }
    return m_Stencil.compute(x);
  }

};

}

#endif // INTERPOLATION_H
