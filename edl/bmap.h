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

#ifndef BMAP_H
#define BMAP_H

#include "edl/edl.h"

namespace EDL_NAMESPACE
{
template<class TValue, class TIndex, int DIM> class BMap;
}

#include "edl/tmdimlist.h"

namespace EDL_NAMESPACE
{

template<class TValue, class TIndex, int DIM>
class BMap
{
public:

  static long int relativeIndex(TMDimList<TValue, TIndex, DIM, BMap>* /** list */,
                                long int i)
  {
    return i;
  }
  
  static long int totalIndex(TMDimList<TValue, TIndex, DIM, BMap> *list,
                             int block, long int i)
  {
    return block * list->maxNumEntries() + i;
  }

  static long int totalIndex(long int max_num_entries, int /* num_blocks */,
                             int block, long int i)
  {
    return block * max_num_entries + i;
  }
};
  
} // namespace

#endif









