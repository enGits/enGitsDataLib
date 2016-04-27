// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2015 enGits GmbH                                         +
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

#ifndef TSDIMVAR_H
#define TSDIMVAR_H

#include "edl/edl.h"

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

  TSDimVar(TMDimList<TValue, TIndex, DIM, MAP> *a_mdim_list, 
		  TMDimIndex<TIndex> an_index);
  TSDimVar() : TMappedVar<TValue, TIndex, DIM, MAP>() {}
  TSDimVar(const TSDimVar<TValue, TIndex, DIM, MAP> &other);
  virtual ~TSDimVar() {}
  virtual void update();

  TValue& operator[](size_t i) const {
#ifdef EDL_DEBUG
    if (!TMappedVar<TValue, TIndex, DIM, MAP>::initialized) {
      std::cerr << "trying to use [] on an uninitialized variable" << std::endl;
      exit(EXIT_FAILURE);
    }
    if ((i < 0) || (i >= TMappedVar<TValue, TIndex, DIM, MAP>::mdim_list->NumEntries())) {
      std::cerr << "index " << i << "out of bounds" << std::endl;
      throw InvalidIndex_error(i);
    }
    if (!m_MDimList->IsActive(i)) {
      std::cerr << "the entry number " << i << "is inactive" << std::endl;
      throw InvalidIndex_error(i);
    }
#endif
    return m_Value[MAP::relativeIndex(TMappedVar<TValue, TIndex, DIM, MAP>::m_MDimList, i)];
  }
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










