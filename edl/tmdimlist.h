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

#ifndef TMDIMLIST_H
#define TMDIMLIST_H

namespace EDL_NAMESPACE
{
template <class TValue, class TIndex, int DIM, class MAP> class TMDimList;
}

#include "edl/tmdimindex.h"
#include "edl/tmappedvar.h"
#include "edl/tlist.h"
#include "edl/tkeylist.h"

namespace EDL_NAMESPACE
{

template<class TValue, class TIndex, int DIM, class MAP>
class TMDimList : public List
{

  friend class TMappedVar<TValue, TIndex, DIM, MAP>;

protected:

  TKeyList<TMDimIndex<TIndex>, TValue*> *m_TotalIndex;

  TIndex          m_EmptyIndexValue;
  TList<TIndex> **m_SubIndex;
  TValue         *m_RootPointer;
  void          **m_VarPointer;
  TValue          m_DefaultValue;

  void AllocNewBlocks();
  virtual void Extend(size_t delta);
  virtual void InitEntry (size_t i);
  TList<TMappedVar<TValue, TIndex, DIM, MAP>*> *m_MappedVars;


public:

  TMDimList(size_t a_max_num_entries, size_t a_delta_entries, TValue a_default_value, TIndex an_empty_index_value);
  TMDimList(List *a_master, TValue a_default_value, TIndex an_empty_index_value);
  virtual ~TMDimList();

  int LengthOfPointerField() const;
  void SetPointerField();
  void CreateTotalIndex(int level, TMDimIndex<TIndex> &I);
  void AddIndex(int i, TIndex new_index);
  TList<TMDimIndex<TIndex> >* GetSubIndices(TMDimIndex<TIndex> I) const;
  virtual void CopyEntry(size_t src, size_t dst);
  TValue* Resolve(TMDimIndex<TIndex> I);
  size_t numBlocks() const { return m_TotalIndex->numEntries(); }
  size_t numSubIndices(int level) const { return m_SubIndex[level]->numEntries(); }
  TIndex subIndex(int level, int i) { size_t si = i; return m_SubIndex[level]->at(si); }
  bool IndexExists(size_t level, TIndex an_index);
};

//
//.. constructor
//
template<class TValue, class TIndex, int DIM, class MAP>
TMDimList<TValue, TIndex, DIM, MAP>::TMDimList(size_t a_max_num_entries, 
                                               size_t a_delta_entries,
                                               TValue a_default_value,
                                               TIndex an_empty_index_value)
  : List(a_max_num_entries, a_delta_entries)
{
  int i;
  m_DefaultValue = a_default_value;
  m_EmptyIndexValue = an_empty_index_value;
  TMDimIndex<TIndex> default_index;
  m_TotalIndex = new TKeyList<TMDimIndex<TIndex>, TValue*>(5, 5, default_index, NULL);
  m_SubIndex = new TList<TIndex>* [DIM];
  for (i = 0; i < DIM; i++) {
    m_SubIndex[i] = new TList<TIndex>(5, 5, m_EmptyIndexValue);
  }
  m_RootPointer = new TValue[numBlocks()];
  m_MappedVars = new TList<TMappedVar<TValue, TIndex, DIM, MAP>*>(10, 10, NULL);
}

template<class TValue, class TIndex, int DIM, class MAP>
TMDimList<TValue, TIndex, DIM, MAP>::TMDimList(List *a_master,
                                               TValue a_default_value,
                                               TIndex an_empty_index_value)
  : List(a_master)
{
  int i;
  m_DefaultValue = a_default_value;
  m_EmptyIndexValue = an_empty_index_value;
  TMDimIndex<TIndex> default_index;
  m_TotalIndex = new TKeyList<TMDimIndex<TIndex>, TValue*>(5, 5, default_index, NULL);
  m_SubIndex = new TList<TIndex>* [DIM];
  for (i = 0; i < DIM; i++) {
    m_SubIndex[i] = new TList<TIndex>(5, 5, m_EmptyIndexValue);
  }
  m_RootPointer = new TValue[numBlocks()];
  m_MappedVars = new TList<TMappedVar<TValue, TIndex, DIM, MAP>*>(10, 10, NULL);
}

template<class TValue, class TIndex, int DIM, class MAP>
TMDimList<TValue, TIndex, DIM, MAP>::~TMDimList()
{
  delete m_TotalIndex;
  for (size_t i = 0; i < DIM; i++) {
    delete m_SubIndex[i];
  }
  delete [] m_SubIndex;
  delete m_MappedVars;
  delete [] m_RootPointer;
}

//
//.. LengthOfPointerField
//
template<class TValue, class TIndex, int DIM, class MAP>
int TMDimList<TValue, TIndex, DIM, MAP>::LengthOfPointerField() const
{
  int i, L, prod;
  L = m_SubIndex[0]->NumEntries();
  prod = L;
  for (i = 1; i < DIM; i++) {
    prod *= m_SubIndex[i]->NumEntries();
    L += prod;
  }
  return L;
}

//
//.. SetPointerField
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::SetPointerField()
{
  void **p = m_VarPointer;
  int prod = 1;
  int i, j;
  for (i = 0; i < DIM-2; i++) {
    prod *= m_SubIndex[i]->NumEntries();
    for (j = 0; j < prod; j++) {
      p[j] = p + prod + j * m_SubIndex[i+1]->NumEntries();
    }
    p += prod;
  }
}

//
//.. AddIndex
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::AddIndex(int i, TIndex new_index)
{
  TMDimIndex<TIndex> I(DIM);
  m_SubIndex[i]->newEntry() = new_index;
  CreateTotalIndex(0, I);
  AllocNewBlocks();
}

//
//.. CreateTotalIndex
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::CreateTotalIndex(int level, TMDimIndex<TIndex> &I)
{
  if (level < DIM) {
    for (size_t i = m_SubIndex[level]->beginIdx(); i < m_SubIndex[level]->endIdx(); i = m_SubIndex[level]->nextIdx(i)) {
      TMDimIndex<TIndex> next_I = I;
      next_I[level] = m_SubIndex[level]->at(i);
      CreateTotalIndex(level + 1, next_I);
    }
  } else {
    try {
      m_TotalIndex->find(I);
    } 
    catch (NotFound_error) {
      m_TotalIndex->newEntry(I) = NULL;
    }
    //if (m_TotalIndex->Find(I) == -1) m_TotalIndex->NewEntry(I) = NULL;
  };
}

//
//.. AllocNewBlocks
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::AllocNewBlocks()
{
  int num_new_blocks = 0;
  TValue *new_m_RootPointer = new TValue[numBlocks() * maxNumEntries()];

  // find out hom many new blocks are needed
  FORALL(i, m_TotalIndex->) {
    if (m_TotalIndex->viAt(i) == NULL) num_new_blocks++;
  }
                                       
  // copy existing data
  //size_t h;
  for (size_t block = 0; block < numBlocks() - num_new_blocks; block++) {
    for (size_t i = 0; i < maxNumEntries(); i++) {
      //h = MAP::TotalIndex(this, block, i);
      new_m_RootPointer[MAP::totalIndex(maxNumEntries(), numBlocks(), block, i)] =
        m_RootPointer[MAP::totalIndex(maxNumEntries(), numBlocks() - num_new_blocks, block, i)];
    }
  }

  // set recently allocated data to the default value
  for (size_t block = numBlocks() - num_new_blocks; block < numBlocks(); block++) {
    for (size_t i = 0; i < maxNumEntries(); i++) {
      new_m_RootPointer[MAP::totalIndex(this, block, i)] = m_DefaultValue;
    }
  }

  // take the new pointer
  delete [] m_RootPointer;
  m_RootPointer = new_m_RootPointer;

  // set the entries in m_TotalIndex to the correct positions in the new field
  FORALL(block, m_TotalIndex->) {
    m_TotalIndex->viAt(block) = m_RootPointer + MAP::totalIndex(this, block, 0);
  }

  // update all connected variables
  FORALL(i, m_MappedVars->) {
    m_MappedVars->at(i)->update();
  }
}

//
//.. Extend
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::Extend(size_t delta)
{
  TValue* new_m_RootPointer = new TValue[numBlocks() * (maxNumEntries() + delta)];

  // find out how many entries need to be copied
  size_t num_copy_entries = maxNumEntries();
  if (delta < 0) num_copy_entries -= delta;

  // copy existing data to the new field
  for(size_t block = m_TotalIndex->beginIdx();
      block < m_TotalIndex->endIdx();
      block = m_TotalIndex->nextIdx(block)) {
    for (size_t i = 0; i < num_copy_entries; i++) {
      new_m_RootPointer[MAP::totalIndex(maxNumEntries() + delta, numBlocks(), block, i)] =
        m_RootPointer[MAP::totalIndex(maxNumEntries(), numBlocks(), block, i)];
    }
  }

  // take the new pointer
  delete [] m_RootPointer;
  m_RootPointer = new_m_RootPointer;

  // set new data to the default value
  for(size_t block = m_TotalIndex->beginIdx();
      block < m_TotalIndex->endIdx();
      block = m_TotalIndex->nextIdx(block)) {
    for (size_t i = num_copy_entries; i < maxNumEntries() + delta; i++) {
      m_RootPointer[MAP::totalIndex(maxNumEntries() + delta, numBlocks(), block, i)]
        = m_DefaultValue;
    }
  }

  // set the entries in m_TotalIndex to the correct positions in the new field
  for(size_t block = m_TotalIndex->beginIdx();
      block < m_TotalIndex->endIdx();
      block = m_TotalIndex->nextIdx(block)) {
    m_TotalIndex->viAt(block) =
      m_RootPointer + MAP::totalIndex(maxNumEntries() + delta, numBlocks(), block, 0);
  }

  // update all connected variables
  for(size_t i = m_MappedVars->beginIdx(); i < m_MappedVars->endIdx(); i = m_MappedVars->nextIdx(i)) {
    m_MappedVars->at(i)->update();
  }
}

//
//.. CopyEntry
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::CopyEntry(size_t src, size_t dst)
{
  List::copyEntry(src, dst);
  FORALL(i, m_TotalIndex->) {
    m_RootPointer[i * maxNumEntries() + dst] = m_RootPointer[i * maxNumEntries() + src];
  }
}

//
//.. InitEntry
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::InitEntry(size_t i)
{
  FORALL(j, m_TotalIndex->) {
    m_TotalIndex->viAt(j)[MAP::relativeIndex(this, i)] = m_DefaultValue;
  }
}

//
//.. GetSubIndices
//
template<class TValue, class TIndex, int DIM, class MAP>
TList<TMDimIndex<TIndex> >* 
TMDimList<TValue, TIndex, DIM, MAP>::GetSubIndices(TMDimIndex<TIndex> I) const
{
  TList<TMDimIndex<TIndex> > *ind;
  if (I.Dim() == DIM) {
    ind = new TList<TMDimIndex<TIndex> >(1, 1);
    ind->NewEntry() = I;
  } else {
    int level = I.Dim();
    int i, j;
    ind = new TList<TMDimIndex<TIndex> >(5, 5);
    TList<TMDimIndex<TIndex> > *sub_ind;    
    FORALL(i, m_SubIndex[level]->) {
      sub_ind = GetSubIndices(I + m_SubIndex[level]->At(i));
      FORALL(j, sub_ind->) ind->NewEntry() = sub_ind->At(j);
      delete sub_ind;
    }
  }
  return ind;
}

//
//.. Resolve
//
template<class TValue, class TIndex, int DIM, class MAP>
TValue* TMDimList<TValue, TIndex, DIM, MAP>::Resolve(TMDimIndex<TIndex> I)
{
  if (I.Dim() != DIM) {
    cerr << "list and index dimension must be identical" <<  endl;
    exit(EXIT_FAILURE);
  } else {
    return m_TotalIndex->at(I);
  }
}

//
//.. IndexExists
//
template<class TValue, class TIndex, int DIM, class MAP>
bool TMDimList<TValue, TIndex, DIM, MAP>::IndexExists(size_t level, TIndex an_index)
{
  bool found = false;
  FORALL(i, m_SubIndex[level]->) {
    if (m_SubIndex[level]->at(i) == an_index) found = true;
  }
  return found;
}

} // namespace

#endif









