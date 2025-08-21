// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef TGRAPH_H
#define TGRAPH_H

#include "edl/edl.h"

namespace EDL_NAMESPACE
{
template <typename T1, typename T2> class TGraph;
}

#include "edl/tsparsetwodimarray.h"

namespace EDL_NAMESPACE
{

template <typename T1, typename T2>
class TGraph : public TSparseTwoDimArray<int>
{

protected: // attributes

  TSparseTwoDimArray<T1> m_LinkPayload;
  TList<T2>              m_Payload;


public: // methods

  /** constructor.
   *  Creates a new TGraph as master and initialises it.
   *  @param mne initial maximum number of entries
   *  @param delta initial increment of entries
   */
  TGraph(size_t mne, size_t delta)
  {
    m_Payload.link(this, "payload");
    m_LinkPayload.linkSparseTwoDimArray(this, "link-payload");
  }

};


} // EDL_NAMESPACE

#endif // TGRAPH_H


