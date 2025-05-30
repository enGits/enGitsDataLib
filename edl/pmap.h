// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef PMAP_H
#define PMAP_H

#include "edl/edl.h"

namespace EDL_NAMESPACE
{
template<class TValue, class TIndex, int DIM> class PMap;
}

#include "edl/tmdimlist.h"

namespace EDL_NAMESPACE
{

template<class TValue, class TIndex, int DIM>
class PMap
{
public:

  ensure_forceinline static long int relativeIndex(TMDimList<TValue, TIndex, DIM, PMap> *list, long int i)
  {
    return list->NumBlocks() * i;
  }

  ensure_forceinline static long int totalIndex(TMDimList<TValue, TIndex, DIM, PMap> *list, int block, long int i)
  {
    return block + list->NumBlocks() * i;
  }

  ensure_forceinline static long int totalIndex(long int max_num_entries, int num_blocks, int block, long int i)
  {
    return block + num_blocks * i;
  }

};

} // namespace

#endif









