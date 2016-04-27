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

#ifndef SPARSETWODIMARRAY_H
#define SPARSETWODIMARRAY_H

#include "edl/edl.h"

namespace EDL_NAMESPACE
{
class SparseTwoDimArray;
}

#include "edl/tlist.h"

namespace EDL_NAMESPACE
{

/** A twodimensional array with a flexible number of entries
 *  in the second dimension.
 *  This data structure uses two Lists to store the two dimensional
 *  information. The first List represents the first dimension and
 *  every entry points to the start of a small block in the second List
 *  where the real data entries are stored. The following sketch tries
 *  to illustrate this:
 *  <pre>
 *  [-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-]
 *   ^        ^           ^        ^        ^     ^     ^        ^                        
 *   |   _____|           |        |        |     |     |        |                            
 *   |  |   ______________|        |        |     |     |        |            
 *   |  |  |   ____________________|        |     |     |        |                       
 *   |  |  |  |   __________________________|     |     |        |                          
 *   |  |  |  |  |   _____________________________|_____|________|                       
 *   |  |  |  |  |  |   __________________________|_____|                                
 *   |  |  |  |  |  |  |   _______________________|                                      
 *   |  |  |  |  |  |  |  |                                                              
 *  [-][-][-][-][-][-][-][-]
 *  </pre>
 *  In SparseTwoDimArray and its child classes the indices will be used
 *  as follows:
 *  <ul>
 *  <li><b>i</b><br>
 *  <i>i</i> is always the index in the first dimension.
 *  In child classes representing matrices this also represents
 *  the row of the matrix. A loop over all entries can be done as
 *  it is known from List:
 *  \code
 *  TSparseTwoDimArray<int> A(10,10);
 *  // ...
 *  // ... array is filled with life
 *  // ...
 *  for(size_t i = A.begin(); i < A.end(); i = A.next(i)) {
 *    // ...
 *    // do something in here
 *    // ...
 *  };
 *  \endcode
 *  Or, if you use the macro version:
 *  \code
 *  TSparseTwoDimArray<int> A(10,10);
 *  // ...
 *  // ... array is filled with life
 *  // ...
 *  FORALL(i, A.) {
 *    // ...
 *    // do something in here
 *    // ...
 *  };
 *  \endcode
 *
 *  </li>
 *  <li><b>j</b><br>
 *  <i>j</i> is only used in child classes, which represent matrices.
 *  <i>j</i> represents the column of a matrix.<br>
 *  </li>
 *  <li><b>k</b><br>
 *  <i>k</i> is the index within a sub-block of the long data list.
 *  A small loop over all entries with a fixed first index of <i>i=5</i>
 *  could be done like this:
 *  \code
 *  TSparseTwoDimArray<real> A(10,10);
 *  // ...
 *  // ... array is filled with life
 *  // ...
 *  for (size_t k = 0; k < Count(5); k++) {
 *    // ...
 *    // do something with the entry (5,k) in here
 *    // ...
 *  };
 *  \endcode<br>
 *  </li>
 *  <li><b>l</b><br>
 *  <i>l</i> is used as an absolute index in the long data list.
 *  </li>
 *  </ul>
 *  
 */
class SparseTwoDimArray : public List
{

private:

  /** The master SparseTwoDimArray for this object.
   */
  SparseTwoDimArray *m_MasterArray;

  /** an array storing the number of non zeros per line.
   */
  TList<size_t> *m_Count;

  /** an array storing the position of the first non zero entry 
   *  in the 'data' field.
   */
  TList<size_t> *m_First;

  /** a dummy placeholder array for the actual data of the entries.
   *  This is not of the same dimension as 'count'!
   */
  List *m_DummyData;


protected:

  /** Get a pointer to the dummy_data list.
   *  @return a pointer to the dummy_data list
   */
  List* dummyData();

  /** Get the index in the data entry for a (i,k) pair.
   *  @param i the index in the first dimension
   *  @param k the index in the second dimension 
   */
  size_t dataIndex(size_t i, size_t k) const;

  /** Create a TList which is linked to the long data List.
   *  @param list this will be a pointer to the created list
   *  @param a_default_value is the default value for the new TList
   */
  template <class TT>
  void createDataList(TList<TT>* &list, TT a_default_value);

  /** Get the position of the first non zero entry.
   *  @param i the index in the first dimension
   *  @return the position of the first entry in the second dimension.
   */
  size_t first(size_t i) const { return (*m_MasterArray->m_First)[i]; }


public:

  /** constructor. 
   *  Creates a new SparseTwoDimArray as master and initialises it.
   *  @param mne initial maximum number of entries
   *  @param delta initial increment of entries 
   */
  SparseTwoDimArray (size_t mne, size_t delta);

  /** constructor. 
   *  Creates a new SparseTwoDimArray which is linked to an existing 
   *  master List and initialises it.
   *  @param a_master the master List to link the new SparseTwoDimArray to
   */
  SparseTwoDimArray (List* a_master);

  /** constructor. 
   *  Creates a new SparseTwoDimArray which is linked to an existing 
   *  master SparseTwoDimArray and initialises it.
   *  @param a_master_array the master SparseTwoDimArray to link the 
   *  new SparseTwoDimArray to
   */
  SparseTwoDimArray (SparseTwoDimArray *a_master_array);
  
  /** destructor.
   */
  virtual ~SparseTwoDimArray ();

  /** Link to an existing master SparseTwoDimArray and initialise it.
   *  @param a_master_array the master SparseTwoDimArray to link the SparseTwoDimArray to
   */
  void linkSparseTwoDimArray(SparseTwoDimArray *a_master_array, std::string link_name);

  /** Get the number of entries in the second dimension.
   *  @param i the index in the first dimension
   *  @return the number of entries in the second dimension
   */
  size_t count(size_t i) const;

  /** Get the total number of entries in a row.
   *  This will only differ from Count(i) in case of symmetric matrix classes.
   *  @param i the row index
   */
  size_t totalCount(size_t i) const;

  /** the main part of the construction process.
   *  This method will be called by all constructors to
   *  setup the structure of the SparseTwoDimArray.
   */
  void setup();

  /** Add an entry in the second dimension.
   *  If the specified first dimension index is not active
   *  an InvalidIndex_error will be thrown.
   *  @param i the index in the first dimension
   *  @return the index in the second dimension.
   */
  size_t addSubEntry(size_t i);

  /** Get the master array of this SparseTwoDimArray.
   *  @return a pointer to the master array of this SparseTwoDimArray
   */
  SparseTwoDimArray* masterArray();

  /** delete a sublist. 
   *  @param i the index in the first dimension
   */
  void delSubList(size_t i);

  /** Resize a sub-array.
   *  This method deletes a subarray and creates a new one with a given length.
   *  @param i the index in the first dimension
   *  @param l the length of the new subarray 
   */
  void resize (size_t i, size_t l);

  /**
   * reset the data array.
   * @param mne the initial maximum number of entries
   * @param delta a delta entries to extend or shrink the list
   */
  virtual void resetData(size_t mne, size_t delta = 0);

  /**
   * reset the data array.
   * @param rel_mne the initial maximum number of entries (multiplied by first dimension)
   * @param delta a delta entries to extend or shrink the list (multiplied by first dimension)
   */
  virtual void relativeResetData(real mne, real delta);

  virtual void deleteData(size_t i);
  virtual void cleanUp();

protected:

  virtual void customReset(real, real) { relativeResetData(2,1); }

};

inline size_t SparseTwoDimArray::count(size_t i) const
{
  if (isActive(i)) {
    return m_MasterArray->m_Count->at(i);
  } else {
    return 0;
  }
}

inline size_t SparseTwoDimArray::totalCount(size_t i) const
{ 
  return count(i);
}

template <class TT>
inline void SparseTwoDimArray::createDataList(TList<TT>* &list, TT a_default_value)
{ 
  list = new TList<TT>(m_DummyData, a_default_value);
}

inline size_t SparseTwoDimArray::dataIndex(size_t i, size_t k) const
{ 
#ifdef EDL_DEBUG
  if (count(i) <= k) {
    std::cerr << "index " << k << "out of bounds" << std::endl;
    throw InvalidIndex_error(k);
  };
#endif
  return first(i) + k;
}

} // namspace

#endif





