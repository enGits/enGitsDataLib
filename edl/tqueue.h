// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2015-2019 enGits GmbH                                    +
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

#ifndef TQUEUE_H
#define TQUEUE_H

#include "edl/tlist.h"

namespace EDL_NAMESPACE
{
template <class T> class TQueue;
template <class T, size_t MNE, size_t DE> class TTQueue;
}

#include <typeinfo>

namespace EDL_NAMESPACE
{

template <typename T>
class TQueue : public TList<int>
{

protected: // attributes

  edl::TList<T> *m_ValueList;
  int            m_Oldest;
  int            m_Newest;
  int            m_Size;


public: // methods

  /**
   * @param a_max_num_entries the initial maximal number of entries
   * @param a_delta_entries the increment for expansion
   * @param a_default_value for new entries
   */
  TQueue(size_t a_max_num_entries , size_t a_delta_entries , T a_default_value = T());

  ~TQueue();

  size_t insert(const T& value);
  T      remove();
  T      value(size_t i) { return m_ValueList->at(i); }
  int    oldestIndex() { return m_Oldest; }
  int    newestIndex() { return m_Newest; }
  T      oldestValue() { return m_ValueList->at(m_Oldest); }
  T      newestValue() { return m_ValueList->at(m_Newest); }
  int    size()        { return m_Size; }

};


template <typename T>
TQueue<T>::TQueue(size_t a_max_num_entries, size_t a_delta_entries, T a_default_value)
  : TList<int> (a_max_num_entries, a_delta_entries, -1)
{
  for (size_t i = 0; i < maxNumEntries(); ++i) m_Value[i] = -1;
  m_ValueList = new TList<T>(this, a_default_value);
  m_Oldest = -1;
  m_Newest = -1;
  m_Size   =  0;
}

template <typename T>
TQueue<T>::~TQueue()
{
  delete m_ValueList;
}

template <typename T>
size_t TQueue<T>::insert(const T& value)
{
  size_t i = addEntry();
  if (m_Newest >= 0) {
    at(m_Newest) = i;
  }
  at(i) = -1;
  m_Newest = i;
  if (m_Oldest < 0) {
    m_Oldest = i;
  }
  m_ValueList->at(i) = value;
  ++m_Size;
  return i;
}

template <typename T>
T TQueue<T>::remove()
{
  T value = m_ValueList->at(m_Oldest);
  int oldest = at(m_Oldest);
  delEntry(m_Oldest);
  m_Oldest = oldest;
  --m_Size;
  return value;
}


template<typename T, size_t MNE, size_t DE>
class TTQueue : public TQueue<T>
{
public:
  TTQueue() : TQueue<T>(MNE, DE, -1) {}
};

}

#endif // TQUEUE_H
