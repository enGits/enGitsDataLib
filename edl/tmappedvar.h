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

#ifndef TMAPPEDVAR_H
#define TMAPPEDVAR_H

#include "edl/edl.h"

namespace EDL_NAMESPACE
{
template<class TValue, class TIndex, int DIM, class MAP> class TMappedVar;
}

#include "edl/tmdimlist.h"

namespace EDL_NAMESPACE
{
template<class TValue, class TIndex, int DIM, class MAP>
class TMappedVar
{
protected:

  TMDimList<TValue, TIndex, DIM, MAP> *m_MDimList;
  TMDimIndex<TIndex> m_Index;
  bool m_Initialised;

public:

  TMappedVar(TMDimList<TValue, TIndex, DIM, MAP> *a_mdim_list, TMDimIndex<TIndex> an_index);
  TMappedVar();
  virtual ~TMappedVar();
  virtual void update() = 0;
  virtual void operator=(const TMappedVar<TValue, TIndex, DIM, MAP> &other);
};


template<class TValue, class TIndex, int DIM, class MAP>
TMappedVar<TValue, TIndex, DIM, MAP>::TMappedVar
(TMDimList<TValue, TIndex, DIM, MAP> *a_mdim_list, TMDimIndex<TIndex> an_index)
{
  m_MDimList = a_mdim_list;
  m_Index = an_index;
  if (m_MDimList != NULL) {
    m_MDimList->m_MappedVars->newEntry() = this;
  }
  m_Initialised = false;
}

template<class TValue, class TIndex, int DIM, class MAP>
TMappedVar<TValue, TIndex, DIM, MAP>::TMappedVar()
{
  m_MDimList = NULL;
  m_Index = TMDimIndex<TIndex>();
  m_Initialised = false;
}

template<class TValue, class TIndex, int DIM, class MAP>
TMappedVar<TValue, TIndex, DIM, MAP>::~TMappedVar()
{
  if (m_MDimList != NULL) {
    try {
      size_t i = m_MDimList->m_MappedVars->findItem(this);
      m_MDimList->m_MappedVars->delEntry(i);
    } catch (NotFound_error) {
      cerr << "error in TMappedVar logic" << endl;
      exit(EXIT_FAILURE);
    }
  }
}

template<class TValue, class TIndex, int DIM, class MAP>
void TMappedVar<TValue, TIndex, DIM, MAP>::operator=
(const TMappedVar<TValue, TIndex, DIM, MAP> &other)
{
  if (m_MDimList != NULL) {
    try {
      size_t i = m_MDimList->m_MappedVars->findItem(this);
      m_MDimList->m_MappedVars->delEntry(i);
    } catch (NotFound_error) {
      cerr << "error in TMappedVar logic" << endl;
      exit(EXIT_FAILURE);
    }
  }
  m_MDimList = other.m_MDimList;
  m_Index = other.m_Index;
  if (m_MDimList != NULL) m_MDimList->m_MappedVars->newEntry() = this;
  update();
}

} // namespace

#endif









