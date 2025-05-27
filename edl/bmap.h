// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
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









