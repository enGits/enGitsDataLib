// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef TKEYLIST_H
#define TKEYLIST_H

#include "edl/edl.h"

namespace EDL_NAMESPACE
{
template <class TKey, class TValue> class TKeyList;
}

#include "edl/tlist.h"

namespace EDL_NAMESPACE
{

/**
 * an assosiative array
 */
template<class TKey, class TValue>
class TKeyList : public TList<TKey>
{
protected:

  /** list for the corresponding values. */
  TList<TValue> *m_Value;


public:

  /** constructor. */
  TKeyList(size_t max, size_t delta, TKey default_key, TValue default_value);

  /** @return the index for a given key. */
  size_t find(TKey key) const;

  /** @return the value for a given key. */
  TValue& at(TKey key);

  TValue& viAt(size_t i) { return m_Value->at(i); }
  TKey& kiAt(size_t i) { return TList<TKey>::at(i); }

  /** enter a new entry. */
  TValue& newEntry (TKey key);

};


//
//.. constructor
//
template<class TKey, class TValue>
TKeyList<TKey, TValue>::TKeyList(size_t max, size_t delta, TKey default_key, TValue default_value) :
  TList<TKey> (max, delta, default_key)
{
  m_Value = new TList<TValue> (this, default_value);
}

//
//.. Find
//
template<class TKey, class TValue>
size_t TKeyList<TKey, TValue>::find(TKey key) const
{
  for(size_t i = this->beginIdx(); i < this->endIdx(); i = this->nextIdx(i)) {
    if (TList<TKey>::at(i) == key) return i;
  }
  throw NotFound_error("TKeyList::Find");
}

//
//.. at
//
template<class TKey, class TValue>
TValue& TKeyList<TKey, TValue>::at(TKey key)
{
  for(size_t i = this->beginIdx(); i < this->endIdx(); i = this->nextIdx(i)) {
    if (TList<TKey>::at(i) == key) return m_Value->at(i);
  }
  //return Value->default_value;
  throw NotFound_error("TKeyList::at");
}

//
//.. NewEntry
//
template<class TKey, class TValue>
TValue& TKeyList<TKey, TValue>::newEntry(TKey key)
{
  size_t i = List::addEntry();
  TList<TKey>::at(i) = key;
  return m_Value->at(i);
}

} // namespace

#endif









