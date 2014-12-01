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









