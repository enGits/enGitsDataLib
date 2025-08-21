// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef TSDIMVAR_H
#define TSDIMVAR_H

#include "edl/edl.h"
#include "edl/edlerror.h"

namespace EDL_NAMESPACE
{
template<class TValue, class TIndex, int DIM, class MAP> class TSDimVar;
}

#include "edl/tmappedvar.h"

namespace EDL_NAMESPACE
{
template<class TValue, class TIndex, int DIM, class MAP>
class TSDimVar : public TMappedVar<TValue, TIndex, DIM, MAP>
{
protected:

  TValue *m_Value;

  using TMappedVar<TValue, TIndex, DIM, MAP>::m_MDimList;


public:

  typedef TValue value_type;

  struct proxy_t
  {
    value_type& value;
    bool validated = false;

    proxy_t(value_type* v) : value(*v) 
    {
      validate();
    }

    void validate() 
    {
      if (isnan(value)) {
        throw EdlError("NaN detected");
      }
      if (isinf(value)) {
        throw EdlError("Inf detected");
      }
      validated = true;
    }

    proxy_t& operator=(const value_type& v) 
    { 
      value = v;
      validate();
      return *this; 
    }

    operator value_type&() 
    {
      return value;
    }

    ~proxy_t() 
    { 
      if (!validated) {
        validate();
      }
    }
  };


  TSDimVar(TMDimList<TValue, TIndex, DIM, MAP> *a_mdim_list, TMDimIndex<TIndex> an_index);
  TSDimVar() : TMappedVar<TValue, TIndex, DIM, MAP>() {}
  TSDimVar(const TSDimVar<TValue, TIndex, DIM, MAP> &other);
  virtual ~TSDimVar() {}
  virtual void update();

#ifdef EDL_DEBUG

  TValue& operator[](size_t i) const {
    //return m_Value[MAP::relativeIndex(TMappedVar<TValue, TIndex, DIM, MAP>::m_MDimList, i)];
    return proxy_t(&m_Value[MAP::relativeIndex(TMappedVar<TValue, TIndex, DIM, MAP>::m_MDimList, i)]);
  }
#else
  TValue& operator[](size_t i) const {
    return m_Value[MAP::relativeIndex(TMappedVar<TValue, TIndex, DIM, MAP>::m_MDimList, i)];
  }
#endif

  virtual void operator=(const TSDimVar<TValue, TIndex, DIM, MAP> &other);
  virtual void print();
};

//
//.. constructor
//
template<class TValue, class TIndex, int DIM, class MAP>
TSDimVar<TValue, TIndex, DIM, MAP>::TSDimVar
(TMDimList<TValue, TIndex, DIM, MAP> *a_mdim_list, TMDimIndex<TIndex> an_index)
  : TMappedVar<TValue, TIndex, DIM, MAP>(a_mdim_list, an_index)
{
  update();
}

template<class TValue, class TIndex, int DIM, class MAP>
TSDimVar<TValue, TIndex, DIM, MAP>::TSDimVar
(const TSDimVar<TValue, TIndex, DIM, MAP> &other)
  : TMappedVar<TValue, TIndex, DIM, MAP>(other.m_MDimList, other.m_Index)
{
  if (other.m_Initialised) update();
}

//
//.. Update
//
template<class TValue, class TIndex, int DIM, class MAP>
void TSDimVar<TValue, TIndex, DIM, MAP>::update()
{
  m_Value = TMappedVar<TValue, TIndex, DIM, MAP>::m_MDimList->resolve(TMappedVar<TValue, TIndex, DIM, MAP>::m_Index);
  TMappedVar<TValue, TIndex, DIM, MAP>::m_Initialised = true;
}

//
//.. operator=
//
template<class TValue, class TIndex, int DIM, class MAP>
void TSDimVar<TValue, TIndex, DIM, MAP>::operator=
(const TSDimVar<TValue, TIndex, DIM, MAP> &other)
{
  TMappedVar<TValue, TIndex, DIM, MAP>::operator=(other);
  update();
}

//
//.. operator<<
//
template<class TValue, class TIndex, int DIM, class MAP>
void TSDimVar<TValue, TIndex, DIM, MAP>::print()
{
  std::cout << TMappedVar<TValue, TIndex, DIM, MAP>::m_Index;
}
  
} // namespace

#endif










