// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef TMDIMVAR_H
#define TMDIMVAR_H

#include "edl/edl.h"
#include "edl/edlerror.h"

namespace EDL_NAMESPACE
{
template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE> class TMDimVar;
}

#include "edl/tmdimindex.h"
#include "edl/tmdimlist.h"
#include "edl/tmappedvar.h"

namespace EDL_NAMESPACE
{
template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
class TMDimVar : public TMappedVar<TValue, TIndex, DIM, MAP>
{
protected:

  MAPVALUE *m_Value;

public:

  TMDimVar(TMDimList<TValue, TIndex, DIM, MAP> *a_mdim_list, TMDimIndex<TIndex> an_index);
  TMDimVar(const TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE> &other);
  TMDimVar();
  virtual ~TMDimVar() { delete [] m_Value; }
  virtual void update();

  const MAPVALUE& operator[](size_t i) const {
    return m_Value[i];
  }

  virtual void operator=(const TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE> &other);
  virtual void print();
};



//
//.. constructor
//
template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE>::TMDimVar
(TMDimList<TValue, TIndex, DIM, MAP> *a_mdim_list, TMDimIndex<TIndex> an_index)
  : TMappedVar<TValue, TIndex, DIM, MAP>(a_mdim_list, an_index)
{
  // check if index is legal
  if (TMappedVar<TValue, TIndex, DIM, MAP>::m_Index.Dim() > DIM) {
    std::cerr << "illegal index dimension for TMDimVar" << std::endl;
    exit(EXIT_FAILURE);
  };
  m_Value = NULL;
  update();
}

template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE>::TMDimVar
(const TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE> &other)
  : TMappedVar<TValue, TIndex, DIM, MAP>(other.m_MDimList, other.m_Index)
{
  m_Value = NULL;
  if (other.m_Initialised) update();
}

template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE>::TMDimVar() 
  : TMappedVar<TValue, TIndex, DIM, MAP>()
{
  m_Value = NULL;
  //??Update();
}

//
//.. Update
//
template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
void TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE>::update()
{
  int i;
  int level = TMappedVar<TValue, TIndex, DIM, MAP>::m_Index.Dim();
  int N = TMappedVar<TValue, TIndex, DIM, MAP>::m_MDimList->numSubIndices(level);
  delete [] m_Value;
  m_Value = new MAPVALUE [N];
  for(i = 0; i < N; i++) {
    m_Value[i] = MAPVALUE(TMappedVar<TValue, TIndex, DIM, MAP>::m_MDimList,
                        TMappedVar<TValue, TIndex, DIM, MAP>::m_Index
                        + TMappedVar<TValue, TIndex, DIM, MAP>::m_MDimList->subIndex(level, i));
  }
  TMappedVar<TValue, TIndex, DIM, MAP>::m_Initialised = true;
}

//
//.. operator=
//
template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
void TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE>::operator=
(const TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE> &other)
{
  TMappedVar<TValue, TIndex, DIM, MAP>::operator=(other);
  update();
}

//
//.. print
//
template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
void TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE>::print()
{
  int level = TMappedVar<TValue, TIndex, DIM, MAP>::m_Index.Dim();
  int N = TMappedVar<TValue, TIndex, DIM, MAP>::m_MDimList->numSubIndices(level);
  std::cout << "[" << TMappedVar<TValue, TIndex, DIM, MAP>::m_Index << "->";
  for(int i = 0; i < N; i++) {
    if (i != 0) std::cout << ",";
    m_Value[i].print();
  }
  std::cout << "]";
}

  
} // namespace

#endif










