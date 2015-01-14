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

#ifndef TSPARSETWODIMARRAY_H
#define TSPARSETWODIMARRAY_H

#include "edl/edl.h"

namespace EDL_NAMESPACE
{
template <class T> class TSparseTwoDimArray;
}

#include "edl/sparsetwodimarray.h"

namespace EDL_NAMESPACE
{

/**
 * A twodimensional array with a flexible number of entries
 * in the second dimension.
 */
template <class T>
class TSparseTwoDimArray : public SparseTwoDimArray
{

private:
  /**
   * an array holding the actual data of the entries.
   * This data is not of the same dimension as 'count'!
   */
  TList<T> *m_Data;

protected:

  /// the default value for "zero" entries.
  T m_DefaultValue;

  /**
   * Access to an entry of the 'data' list.
   * @param l the index of the entry
   * @return the value of the entry
   */
  T& DataEntry(size_t l) { return (*m_Data)[l]; }


public:

  /**
   * constructor.
   * Creates a new SparseTwoDimArray as master and initialises it.
   * @param mne initial maximum number of entries
   * @param delta initial increment of entries
   */
  TSparseTwoDimArray (size_t mne, size_t delta);

  /**
   * constructor.
   * Creates a new SparseTwoDimArray as master and initialises it.
   * @param mne initial maximum number of entries
   * @param delta initial increment of entries
   * @param a_default_value a default value for this array
   */
  TSparseTwoDimArray (size_t mne, size_t delta, T a_default_value);

  /**
   * constructor.
   * Creates a new SparseTwoDimArray which is linked to an existing
   * master List and initialises it.
   * @param a_master the master List to link the new SparseTwoDimArray to
   */
  TSparseTwoDimArray (List* a_master);

  /**
   * constructor.
   * Creates a new SparseTwoDimArray which is linked to an existing
   * master List and initialises it.
   * @param a_master the master List to link the new SparseTwoDimArray to
   * @param a_default_value a default value for this array
   */
  TSparseTwoDimArray (List* a_master, T a_default_value);

  /**
   * constructor.
   * Creates a new SparseTwoDimArray which is linked to an existing
   * master SparseTwoDimArray and initialises it.
   * @param a_master_array the master SparseTwoDimArray to link the
   * new SparseTwoDimArray to
   */
  TSparseTwoDimArray (SparseTwoDimArray *a_master_array);

  /**
   * constructor.
   * Creates a new SparseTwoDimArray which is linked to an existing
   * master SparseTwoDimArray and initialises it.
   * @param a_master_array the master SparseTwoDimArray to link the
   * new SparseTwoDimArray to
   * @param a_default_value a default value for this array
   */
  TSparseTwoDimArray (SparseTwoDimArray *a_master_array, T a_default_value);

  /// destructor.
  virtual ~TSparseTwoDimArray ();

  /**
   * Access to an entry.
   * This method provides access to a specific entry of a subarray.
   * @param i the index in the first dimension
   * @param k the index in the second dimension
   * @return a reference to the entry
   */
  T      & at(size_t i, size_t k)       { return (*m_Data)[dataIndex(i,k)]; }
  T const& at(size_t i, size_t k) const { return (*m_Data)[dataIndex(i,k)]; }

  /**
   * Create a new sub entry and return a reference to it.
   * If the specified first dimension index is not active
   * an InvalidIndex_error will be thrown.
   * @param i the index in the first dimension
   * @return a reference to the new entry
   */
  T& newSubEntry(size_t i);

  /// Print matrix to a stream in skyline format.
  void printSkyline(ostream &s, uint_t width = 0);

  /**
   * Get the default value of this array.
   * @return the default value
   */
  T defaultValue() { return m_DefaultValue; }

  /**
   * find an entry in a sub-list.
   * This method will throw a NotFound_error if it does
   * not find the entry.
   * @param i the first index
   * @param a_T the value to look for
   * @return the sub-list index (k)
   */
  size_t findSubItem(size_t i, T a_T);

};


template <class T>
TSparseTwoDimArray<T>::TSparseTwoDimArray(size_t mne, size_t delta)
  : SparseTwoDimArray(mne, delta)
{
  m_Data = new TList<T>(dummyData());
}

template <class T>
TSparseTwoDimArray<T>::TSparseTwoDimArray(size_t mne, size_t delta, T a_default_value)
  : SparseTwoDimArray(mne, delta)
{
  m_DefaultValue = a_default_value;
  m_Data = new TList<T>(dummyData(), m_DefaultValue);
}

template <class T>
TSparseTwoDimArray<T>::TSparseTwoDimArray(List *a_master)
  : SparseTwoDimArray(a_master)
{
  m_Data = new TList<T>(dummyData());
}

template <class T>
TSparseTwoDimArray<T>::TSparseTwoDimArray(List *a_master, T a_default_value)
  : SparseTwoDimArray(a_master)
{
  m_DefaultValue = a_default_value;
  m_Data = new TList<T>(dummyData(), m_DefaultValue);
}

template <class T>
TSparseTwoDimArray<T>::TSparseTwoDimArray(SparseTwoDimArray *a_master_array)
  : SparseTwoDimArray(a_master_array)
{
  m_Data = new TList<T>(dummyData());
}

template <class T>
TSparseTwoDimArray<T>::TSparseTwoDimArray(SparseTwoDimArray *a_master_array, T a_default_value)
  : SparseTwoDimArray(a_master_array)
{
  m_DefaultValue = a_default_value;
  m_Data = new TList<T>(dummyData(), m_DefaultValue);
}

template <class T>
TSparseTwoDimArray<T>::~TSparseTwoDimArray()
{
}

template <class T>
void TSparseTwoDimArray<T>::printSkyline(ostream &s, uint_t width)
{
  FORALL(i, this->) {
    if (width) s.width(4);
    s << i;
    s << "->";
    for (size_t j = 0; j < count(i); j++) {
      if (width) s.width(width);
      s << (*m_Data)[first(i) + j] << ' ';
    }
    s << endl;
  }
}

template <class T>
T& TSparseTwoDimArray<T>::newSubEntry(size_t i)
{
  size_t k = addSubEntry(i);
  return at(i,k);
}

template <class T>
size_t TSparseTwoDimArray<T>::findSubItem(size_t i, T a_T)
{
  size_t k = 0;
  while (k < count(i)) {
    if (at(i, k) == a_T) break;
    k++;
  }
  if (at(i, k) == a_T) return k;
  string msg = "TSparseTwoDimArray::FindSubItem\n";
  msg += typeid(*this).name();
  throw NotFound_error(msg);
}

} // namespace

#endif
