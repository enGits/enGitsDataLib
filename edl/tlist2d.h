// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef TLIST2D_H
#define TLIST2D_H

#include "edl/edl.h"
namespace EDL_NAMESPACE
{
template <class T> class TList2D;
}

#include "edl/list.h"

namespace EDL_NAMESPACE
{

/**
   a twodimensional TList.
*/
template <class T>
class TList2D : public List
{

protected:

  int m_Dim;

  /** array of the stored values */
  T **m_Value;

  virtual void extend (long int);
  virtual void copyEntry (long int source_point, long int dest_point);
  virtual void initEntry (long int i);


public:

  T m_DefaultValue; ///< the default value for new entries

  /** constructor.
    @param a_max_num_entries the initial maximal number of entries
    @param a_delta_entries the initial increment
    @param a_dim the 2nd dimension */
  TList2D(long int a_max_num_entries, long int a_delta_entries, int a_dim);

  /** constructor.
    @param a_max_num_entries the initial maximal number of entries
    @param a_delta_entries the initial increment
    @param a_dim the 2nd dimension
    @param a_default_value the initial value for new entries */
  TList2D(long int a_max_num_entries, long int a_delta_entries,
     int a_dim, T a_default_value);

  /** constructor.
    @param a_dim the 2nd dimension
    @param a_master a master for this List */
  TList2D(int a_dim, List *a_master, std::string link_name = "__none");

  /** constructor.
    @param a_dim the 2nd dimension
    @param a_master a master for this List.
    @param a_default_value the initial value for new entries */
  TList2D(int a_dim, List *a_master, T a_default_value);

  /** destructor. */
  ~TList2D();

  /** The Value of an entry.
    @param i the 1st index of the entry (the shorter one)
    @param j the 2nd index of the entry
    @return the value of this entry */
  T& at (int i, long int j) { return m_Value[i][j]; }

  /** Get a pointer to the data field */
  T** valuePointer() { return m_Value; }

};

//
//.. constructor
//
template <class T>
TList2D<T>::TList2D (long int a_max_num_entries, long int a_delta_entries, int a_dim)
  : List (a_max_num_entries, a_delta_entries)
{
  int i;
  m_Dim = a_dim;
  m_Value = new T* [m_Dim];
  m_Value[0] = new T [maxNumEntries() * m_Dim];
  for (i = 1; i < m_Dim; i++) m_Value[i] = m_Value[i-1] + maxNumEntries();
}

template <class T>
TList2D<T>::TList2D (long int a_max_num_entries, long int a_delta_entries,
                     int a_dim, T a_default_value)
  : List (a_max_num_entries, a_delta_entries)
{
  int i;
  m_Dim = a_dim;
  m_Value = new T* [m_Dim];
  m_Value[0] = new T [maxNumEntries() * m_Dim];
  for (i = 1; i < m_Dim; i++) m_Value[i] = m_Value[i-1] + maxNumEntries();
  m_DefaultValue = a_default_value;
}

template <class T>
TList2D<T>::TList2D (int a_dim, List *a_master, std::string link_name) : List (a_master, link_name)
{
  int i;
  m_Dim = a_dim;
  m_Value = new T* [m_Dim];
  m_Value[0] = new T [maxNumEntries() * m_Dim];
  for (i = 1; i < m_Dim; i++) m_Value[i] = m_Value[i-1] + maxNumEntries();
}

template <class T>
TList2D<T>::TList2D (int a_dim, List *a_master, T a_default_value) : List (a_master)
{
  int i;
  m_Dim = a_dim;
  m_Value = new T* [m_Dim];
  m_Value[0] = new T [maxNumEntries() * m_Dim];
  for (i = 1; i < m_Dim; i++) m_Value[i] = m_Value[i-1] + maxNumEntries();
  m_DefaultValue = a_default_value;
}

//
//.. destructor
//
template <class T>
TList2D<T>::~TList2D ()
{
  delete [] m_Value[0];
  delete [] m_Value;
}

//
//.. Extend
//
template <class T>
void TList2D<T>::extend(long int delta)
{
  int i;
  long int i_entry;
  T **new_value;
  new_value = new T* [m_Dim];
  new_value[0] = new T [(maxNumEntries() + delta) * m_Dim];
  for (i = 1; i < m_Dim; i++) new_value[i] = new_value[i-1] + maxNumEntries() + delta;
  for (i = 0; i < m_Dim; i++) {
    for (i_entry = 0; i_entry < maxNumEntries(); i_entry++) {
      new_value[i][i_entry] = m_Value[i][i_entry];
    }
  }
  delete [] m_Value[0];
  delete [] m_Value;
  m_Value = new_value;
}

//
//.. CopyEntry
//
template <class T>
void TList2D<T>::copyEntry(long int src, long int dst)
{
  int i;
  for (i = 0; i < m_Dim; i++) at(i, dst) = at(i, src);
}

//
//.. InitEntry
//
template <class T>
void TList2D<T>::initEntry(long int j)
{
  int i;
  for (i = 0; i < m_Dim; i++) at(i, j) = m_DefaultValue;
}

} // namespace

#endif
