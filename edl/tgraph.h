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


