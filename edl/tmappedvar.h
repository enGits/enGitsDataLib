// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
  void disable() { m_MDimList = NULL; }

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
      std::cerr << "error in TMappedVar logic" << std::endl;
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
      std::cerr << "error in TMappedVar logic" << std::endl;
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









