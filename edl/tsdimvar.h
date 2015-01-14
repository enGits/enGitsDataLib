// !!
// This is a part of MOUSE, a library for PDE's on unstructured grids
// Copyright (C) 1999 Oliver Gloth <oliver@vug.uni-duisburg.de>
// Institut fuer Verbrennung und Gasdynamik (Universitaet Duisburg)
// institute for combustion and gas dynamics (university of Duisburg)
// Thursday, 1 April, 1999 Duisburg, Germany
//
// please see http://www.vug.uni-duisburg.de/MOUSE for more information
// please send any questions or suggestions to mouse@www.vug.uni-duisburg.de
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software 
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// !!

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
      cerr << "trying to use [] on an uninitialized variable" << endl;
      exit(EXIT_FAILURE);
    }
    if ((i < 0) || (i >= TMappedVar<TValue, TIndex, DIM, MAP>::mdim_list->NumEntries())) {
      cerr << "index " << i << "out of bounds" << endl;
      throw InvalidIndex_error(i);
    }
    if (!m_MDimList->IsActive(i)) {
      cerr << "the entry number " << i << "is inactive" << endl;
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
  m_Value = TMappedVar<TValue, TIndex, DIM, MAP>::m_MDimList->Resolve(TMappedVar<TValue, TIndex, DIM, MAP>::m_Index);
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
  cout << TMappedVar<TValue, TIndex, DIM, MAP>::m_Index;
}
  
} // namespace

#endif










