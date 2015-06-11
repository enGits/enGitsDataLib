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

  static long int relativeIndex(TMDimList<TValue, TIndex, DIM, BMap>*, long int i)
  {
    return i;
  }
  
  static long int totalIndex(TMDimList<TValue, TIndex, DIM, BMap> *list, int block, long int i)
  {
    return block * list->maxNumEntries() + i;
  }

  static long int totalIndex(long int max_num_entries, int, int block, long int i)
  {
    return block * max_num_entries + i;
  }
};
  
} // namespace

#endif









