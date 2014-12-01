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

#ifndef TMDIMVAR_H
#define TMDIMVAR_H

#include "edl/edl.h"

namespace EDL_NAMESPACE
{
template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE> class TMDimVar;
}

#include "edl/tmddimindex.h"
#include "edl/tmdimlist.h"
#include "edl/tmappedvar.h"

namespace EDL_NAMESPACE
{
template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
class TMDimVar : public TMappedVar<TValue, TIndex, DIM, MAP>
{
protected:

  MAPVALUE *value;

public:

  TMDimVar(TMDimList<TValue, TIndex, DIM, MAP> *a_mdim_list, 
		  TMDimIndex<TIndex> an_index);
  TMDimVar(const TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE> &other);
  TMDimVar();
  virtual ~TMDimVar() { delete [] value; };
  virtual void Update();
  const MAPVALUE& operator[](size_t i) const { 
#ifdef MDEBUG
    if (!TMappedVar<TValue, TIndex, DIM, MAP>::initialized) {
      cerr << "trying to use [] on an uninitialized variable" << endl;
      throw InvalidIndex_error(i);
    };      
    //    if ((i < 0) || (i >= TMappedVar<TValue, TIndex, DIM, MAP>::mdim_list->NumSubIndices(TMappedVar<TValue, TIndex, DIM, MAP>::index.Dim()))) {
    if (i >= TMappedVar<TValue, TIndex, DIM, MAP>::mdim_list->NumSubIndices(TMappedVar<TValue, TIndex, DIM, MAP>::index.Dim())) {
      cerr << "index " << i << "out of bounds" << endl;
      throw InvalidIndex_error(i);
    };
#endif
    return value[i]; 
  };
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
  if (TMappedVar<TValue, TIndex, DIM, MAP>::index.Dim() > DIM) {
    cerr << "illegal index dimension for TMDimVar" << endl;
    exit(EXIT_FAILURE);
  };
  value = NULL;
  Update();
}

template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE>::TMDimVar
(const TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE> &other)
  : TMappedVar<TValue, TIndex, DIM, MAP>(other.mdim_list, other.index)
{
  value = NULL;
  if (other.initialized) Update();
}

template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE>::TMDimVar() 
  : TMappedVar<TValue, TIndex, DIM, MAP>()
{
  value = NULL;
  //??Update();
}

//
//.. Update
//
template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
void TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE>::Update()
{
  int i;
  int level = TMappedVar<TValue, TIndex, DIM, MAP>::index.Dim();
  int N = TMappedVar<TValue, TIndex, DIM, MAP>::mdim_list->NumSubIndices(level);
  delete [] value;
  value = new MAPVALUE [N];
  for(i = 0; i < N; i++) {
    value[i] = MAPVALUE(TMappedVar<TValue, TIndex, DIM, MAP>::mdim_list, 
                        TMappedVar<TValue, TIndex, DIM, MAP>::index 
                        + TMappedVar<TValue, TIndex, DIM, MAP>::mdim_list->SubIndex(level, i));
  };
  TMappedVar<TValue, TIndex, DIM, MAP>::initialized = true;
}

//
//.. operator=
//
template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
void TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE>::operator=
(const TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE> &other)
{
  TMappedVar<TValue, TIndex, DIM, MAP>::operator=(other);
  Update();
}

//
//.. print
//
template<class TValue, class TIndex, int DIM, class MAP, class MAPVALUE>
void TMDimVar<TValue, TIndex, DIM, MAP, MAPVALUE>::print()
{
  int level = TMappedVar<TValue, TIndex, DIM, MAP>::index.Dim();
  int N = TMappedVar<TValue, TIndex, DIM, MAP>::mdim_list->NumSubIndices(level);
  cout << "[" << TMappedVar<TValue, TIndex, DIM, MAP>::index << "->";
  for(int i = 0; i < N; i++) {
    if (i != 0) cout << ",";
    value[i].print();
  };
  cout << "]";
}

  
} // namespace

#endif










