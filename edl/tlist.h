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

#ifndef TLIST_H
#define TLIST_H

#include <cstddef>
#include <typeinfo>

#include "edl/list.h"

namespace EDL_NAMESPACE
{
template <class T> class TList;
template <class T, size_t MNE, size_t DE> class TTList;
}

#include <typeinfo>

namespace EDL_NAMESPACE
{

/**
 * a self-allocating one-dimensional array.
 */
template <class T>
class TList : public List 
{
private:

  /// a flag which is true if Has or FindItem has already been successfully called
  bool m_FoundLastTime;

  /// the index of the last item found
  size_t m_LastFound;

protected:

  /** array of the stored values */
  T *m_Value;
  
  /** extends the array value.
      @param delta the increment */ 
  virtual void extend (size_t delta);

  /** initialize an entry with the default value 
      @param i the entry to initialize */ 
  virtual void initEntry (size_t i)
  {
    m_Value[i] = m_DefaultValue;
  }
  
  virtual void copyEntry (size_t source_point, size_t dest_point)
  {
    List::copyEntry(source_point, dest_point);
    m_Value[dest_point] = m_Value[source_point];
  }

  /// the default value for new entries
  T m_DefaultValue;

  virtual size_t dataLength() { return sizeof(T); }

  // QT needs to be removed
  // virtual QByteArray partialBuffer(size_t);
  // virtual void       fromPartialBuffer(size_t, QByteArray);


public:
    
  /**
   * default constructor.
   * @param a_default_value for new entries 
   */
  TList(T a_default_value = T());

  /** 
   * @param a_max_num_entries the initial maximal number of entries
   * @param a_delta_entries the increment for expansion
   * @param a_default_value for new entries 
   */
  TList(size_t a_max_num_entries , size_t a_delta_entries , T a_default_value = T());
    
  /**
   * @param a_master a master for this List.
   * @param a_default_value the initial value for new entries 
   */
  TList(List *a_master, T a_default_value = T(), std::string link_name = "__none");

  TList(const TList<T> &other, const T a_default_value = T());
  
  ~TList ();

  /**
   * @brief Set the default value of this container.
   * @param default_value the new default value
   */
  void setDefaultValue(T default_value) { m_DefaultValue = default_value; }
  
  /**
   * Initialize all active entries with a value.
   * @param v the value for the initilaziation 
   */
  void initAllEntries (T v);

  /**
   * Allocate a new entry.
   * @return a reference to the new entry 
   */
  T& newEntry ()
  {
    size_t i = addEntry ();
    return m_Value[i];
  }

  /**
   * the Value of an entry.
   * @param i the index of the entry
   * @return a reference to the value of this entry 
   */
  T& at (size_t i)
  {
#ifdef EDL_DEBUG
    if ((i < 0) || (i >= endIdx())) {
      std::cerr << "index " << i << "out of bounds" << std::endl;
      throw InvalidIndex_error(i);
    }
    if (!isActive(i)) {
      std::cerr << "the entry number " << i << "is inactive" << std::endl;
      throw InvalidIndex_error(i);
    }
#endif
    return m_Value[i];
  }

  /** The Value of an entry.
    @param i the index of the entry
    @return the value of this entry */
  T const& at (size_t i) const
  {
#ifdef EDL_DEBUG
    if ((i < 0) || (i >= endIdx())) {
      std::cerr << "index " << i << "out of bounds" << std::endl;
      throw InvalidIndex_error(i);
    };
    if (!isActive(i)) {
      std::cerr << "the entry number " << i << "is inactive" << std::endl;
      throw InvalidIndex_error(i);
    };
#endif
    return m_Value[i];
  }

  /** The Value of an entry.
      @param i the index of the entry
      @return a reference to the value of this entry */
  T& operator[] (size_t i) { return at(i); }

  /** The Value of an entry.
      @param i the index of the entry
      @return the value of this entry */
  const T operator[] (size_t i) const { return at(i); }

  /** Find an entry.
    @param v the value to search for
    @return the index of the entry */
  size_t findItem (T v);

  /** Check if an entry is in the list */
  bool has(T v);

  /** Get a pointer to the data field */
  T* data() { return m_Value; }

  /** exchange two entries
   *  @param i1 the first entry
   *  @param i2 the second entry
   */ 
  void swap(size_t i1, size_t i2)
  {
    T h = at(i1);
    at(i1) = at(i2);
    at(i2) = h;
  }

  /** set all entries to a constant value
   *  @param v the value to set
   */ 
  void setAll(T v) { FORALL(i,this->) { at(i) = v; }; }

  /** */ void operator=(const TList<T> &other);

  /** return the first active entry
   *  @return a reference to the first active entry.
   */
  T& firstEntry() { return at(beginIdx()); }

  /** return the first active entry
   *  @return a reference to the last active entry.
   */
  T& lastEntry() { return at(lastIdx()); }

  /** Output operator.
   */ 
  template <class aT> 
  friend std::ostream& operator<<(std::ostream &s, const TList<aT> &tlist);


  // =====================================
  // STL style iterator and related things
  // =====================================

  typedef T value_type;

  class iterator {

    TList<T> *tlist;
    size_t i_current;

  public:
    iterator(TList<T> *a_tlist, size_t i) : tlist(a_tlist), i_current(i) {}
    iterator() : tlist(NULL), i_current(0) {}
    bool operator==(const iterator &iter) const;
    iterator& operator++();
    iterator operator++(int);
    T& operator*();
  };

  class const_iterator {

    const TList<T> *tlist;
    size_t i_current;

  public:
    const_iterator(const TList<T> *a_tlist, size_t i) : tlist(a_tlist), i_current(i) {}
    const_iterator() : tlist(NULL), i_current(0) {}
    bool operator==(const const_iterator &iter) const;
    const_iterator& operator++();
    const_iterator operator++(int);
    T operator*();
  };

  iterator begin() { return iterator(this, beginIdx()); }
  iterator end() { return iterator(this, endIdx()); }
  const_iterator begin() const { return const_iterator(this, beginIdx()); }
  const_iterator end() const { return const_iterator(this, endIdx()); }



  /**
   * sort the TList.
   * This only works for types <it>T</it> which have
   * the operator '>' and '='. The smallest item will be
   * the first after this sorting.
   */
  void sortUp();

  /**
   * sort the TList.
   * This only works for types <it>T</it> which have
   * the operator '<' and '='. The smallest item will be
   * the last after this sorting.
   */
  void sortDown();
};





// TList<T>::iterator
// ==================
template<class T>
inline typename TList<T>::iterator& TList<T>::iterator::operator++()
{ 
  i_current = tlist->next_idx(i_current);
  return *this;
}

template<class T>
inline typename TList<T>::iterator TList<T>::iterator::operator++(int)
{
  size_t old_i_current;
  i_current = tlist->next_idx(i_current);
  return iterator(tlist, old_i_current); 
}

template<class T>
inline bool TList<T>::iterator::operator==(const iterator &iter) const
{
  return ((iter.tlist == tlist) && (iter.i_current == i_current));
}

template<class T>
inline T& TList<T>::iterator::operator*()
{
  return (*tlist)[i_current];
}



// TList<T>::const_iterator
// ========================
template<class T>
inline typename TList<T>::const_iterator& TList<T>::const_iterator::operator++()
{ 
  i_current = tlist->next_idx(i_current);
  return const_iterator(tlist, i_current); 
};

template<class T>
inline typename TList<T>::const_iterator TList<T>::const_iterator::operator++(int)
{
  size_t old_i_current;
  i_current = tlist->next_idx(i_current);
  return const_iterator(tlist, old_i_current); 
}

template<class T>
inline bool TList<T>::const_iterator::operator==(const const_iterator &iter) const
{
  return ((iter->tlist == tlist) && (iter.i_current == i_current));
}

template<class T>
inline T TList<T>::const_iterator::operator*()
{
  return (*tlist)[i_current];
}





template<class T, size_t MNE, size_t DE>
class TTList : public TList<T>
{
public:
  TTList() : TList<T>(MNE, DE) {}
};

template <class T>
TList<T>::TList (T a_default_value) : List()
{
  m_Value = NULL;
  m_DefaultValue = a_default_value;
  m_FoundLastTime = false;
}

template <class T>
TList<T>::TList (size_t a_max_num_entries, size_t a_delta_entries, T a_default_value) 
  : List (a_max_num_entries, a_delta_entries)
{
  size_t mne = maxNumEntries();
  m_Value = new T [mne];
  m_DefaultValue = a_default_value;
  m_FoundLastTime = false;
}

template <class T>
TList<T>::TList (List *a_master, T a_default_value, std::string link_name) : List (a_master, link_name)
{
  size_t mne = maxNumEntries();
  m_Value = new T [mne];
  m_DefaultValue = a_default_value;
  initAllEntries(m_DefaultValue);
  m_FoundLastTime = false;
}

template <class T>
TList<T>::TList(const TList<T> &other, const T a_default_value) : List()
{
  m_Value = NULL;
  m_DefaultValue = a_default_value;
  m_FoundLastTime = false;
  operator=(other);
}

template <class T>
void TList<T>::operator=(const TList<T> &other)
{
  List::operator=(other);
  for(size_t i = beginIdx(); i < endIdx(); i = nextIdx(i)) {
    m_Value[i] = other[i];
  }
}

template <class T>
TList<T>::~TList ()
{
  if (m_Value) delete [] m_Value;
} 

//
//.. InitAllEntries
//
template <class T>
void TList<T>::initAllEntries (T a_value) {
  size_t i;
  for (i = 0; i < maxNumEntries(); i++) m_Value[i] = a_value;
}

//
//.. Extend
//
template <class T>
void TList<T>::extend (size_t delta) {
  T *new_value;
  size_t i_entry;
  new_value = new T [maxNumEntries() + delta];
  for (i_entry = 0; i_entry < maxNumEntries(); i_entry++) new_value[i_entry] = m_Value[i_entry];
  delete [] m_Value;
  m_Value = new_value;
}

//
//.. FindItem
//
template <class T>
size_t TList<T>::findItem (T item)
{
  if (m_FoundLastTime) {
    if (isActive(m_LastFound)) {
      if (m_Value[m_LastFound] == item) return m_LastFound;
    }
  }
  for(size_t i = beginIdx(); i < endIdx(); i = nextIdx(i)) {
    if (m_Value[i] == item) {
      m_FoundLastTime = true;
      m_LastFound = i;
      return i;
    }
  }
  std::string msg = "TList::FindItem\n";
  msg += typeid(*this).name();
  throw NotFound_error(msg);
}

//
//.. IsIn
//
template <class T>
bool TList<T>::has(T item)
{
  if (m_FoundLastTime) {
    if (isActive(m_LastFound)) {
      if (m_Value[m_LastFound] == item) return true;
    }
  }
  for(size_t i = beginIdx(); i < endIdx(); i = nextIdx(i)) {
    if (m_Value[i] == item) {
      m_FoundLastTime = true;
      m_LastFound = i;
      return true;
    }
  }
  return false;
}

template <class T>
std::ostream& operator<<(std::ostream &s, const TList<T> &tlist)
{
  s << '[';
  FORALL(i, tlist.) {
    s << "(" << i << "," << tlist[i] << ")";
    if (i != tlist.last_idx()) s << ',';
  }
  s << ']';
  return s;
}

template <class T>
void TList<T>::sortUp()
{
  size_t i_start = beginIdx();
  while (i_start < endIdx()) {
    for (size_t i = i_start; i < lastIdx(); i = nextIdx(i)) {
      if (at(i) > at(nextIdx(i))) {
        swap(i, nextIdx(i));
      }
    }
    i_start = nextIdx(i_start);
  }
}

template <class T>
void TList<T>::sortDown()
{
  size_t i_start = beginIdx();
  while (i_start < endIdx()) {
    for (size_t i = i_start; i < lastIdx(); i = nextIdx(i)) {
      if (at(i) < at(nextIdx(i))) {
        swap(i, nextIdx(i));
      }
    }
    i_start = nextIdx(i_start);
  }
}

// template <class T>
// QByteArray TList<T>::partialBuffer(size_t i)
// {
//   return QByteArray((char*) (m_Value + i), sizeof(T));
// }

// template <class T>
// void TList<T>::fromPartialBuffer(size_t i, QByteArray buffer)
// {
//   if (buffer.length() != sizeof(T)) {
//     std::cerr << "data corrupted in buffer" << std::endl;
//     exit(EXIT_FAILURE);
//   }
//   m_Value[i] = *((T*) buffer.data());
// }

} //namespace

#endif

