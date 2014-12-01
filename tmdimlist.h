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
  TIndex empty_index_value;
  TKeyList<TMDimIndex<TIndex>, TValue*> *total_index;
  TList<TIndex> **sub_index;
  TValue *root_pointer;
  void **var_pointer;
  void AllocNewBlocks();
  virtual void Extend(size_t delta);
  virtual void InitEntry (size_t i);
  TValue default_value;
  TList<TMappedVar<TValue, TIndex, DIM, MAP>*> *mapped_vars;

public:
  TMDimList(size_t a_max_num_entries, size_t a_delta_entries,
		   TValue a_default_value, TIndex an_empty_index_value);
  TMDimList(List *a_master, TValue a_default_value, 
		   TIndex an_empty_index_value);
  virtual ~TMDimList();

  int LengthOfPointerField() const;
  void SetPointerField();
  void CreateTotalIndex(int level, TMDimIndex<TIndex> &I);
  void AddIndex(int i, TIndex new_index);
  TList<TMDimIndex<TIndex> >* GetSubIndices(TMDimIndex<TIndex> I) const;
  virtual void CopyEntry(size_t src, size_t dst);
  TValue* Resolve(TMDimIndex<TIndex> I);
  size_t NumBlocks() const { return total_index->NumEntries(); };
  size_t NumSubIndices(int level) const { return sub_index[level]->NumEntries(); };
  TIndex SubIndex(int level, int i) { size_t si = i; return sub_index[level]->At(si); };
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
  default_value = a_default_value;
  empty_index_value = an_empty_index_value;
  TMDimIndex<TIndex> default_index;
  total_index = new TKeyList<TMDimIndex<TIndex>, TValue*>(5, 5, default_index, NULL);
  sub_index = new TList<TIndex>* [DIM];
  for (i = 0; i < DIM; i++) {
    sub_index[i] = new TList<TIndex>(5, 5, empty_index_value);
  };
  root_pointer = new TValue[NumBlocks()];
  mapped_vars = new TList<TMappedVar<TValue, TIndex, DIM, MAP>*>(10, 10, NULL);
}

template<class TValue, class TIndex, int DIM, class MAP>
TMDimList<TValue, TIndex, DIM, MAP>::TMDimList(List *a_master,
                                               TValue a_default_value,
                                               TIndex an_empty_index_value)
  : List(a_master)
{
  int i;
  default_value = a_default_value;
  empty_index_value = an_empty_index_value;
  TMDimIndex<TIndex> default_index;
  total_index = new TKeyList<TMDimIndex<TIndex>, TValue*>(5, 5, default_index, NULL);
  sub_index = new TList<TIndex>* [DIM];
  for (i = 0; i < DIM; i++) {
    sub_index[i] = new TList<TIndex>(5, 5, empty_index_value);
  };
  root_pointer = new TValue[NumBlocks()];
  mapped_vars = new TList<TMappedVar<TValue, TIndex, DIM, MAP>*>(10, 10, NULL);
}

template<class TValue, class TIndex, int DIM, class MAP>
TMDimList<TValue, TIndex, DIM, MAP>::~TMDimList()
{
  delete total_index;
  for (size_t i = 0; i < DIM; i++) {
    delete sub_index[i];
  };
  delete [] sub_index;
  delete mapped_vars;
  delete [] root_pointer;
}

//
//.. LengthOfPointerField
//
template<class TValue, class TIndex, int DIM, class MAP>
int TMDimList<TValue, TIndex, DIM, MAP>::LengthOfPointerField() const
{
  int i, L, prod;
  L = sub_index[0]->NumEntries();
  prod = L;
  for (i = 1; i < DIM; i++) {
    prod *= sub_index[i]->NumEntries();
    L += prod;
  };
  return L;
}

//
//.. SetPointerField
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::SetPointerField()
{
  void **p = var_pointer;
  int prod = 1;
  int i, j;
  for (i = 0; i < DIM-2; i++) {
    prod *= sub_index[i]->NumEntries();
    for (j = 0; j < prod; j++) {
      p[j] = p + prod + j * sub_index[i+1]->NumEntries();
    };
    p += prod;
  };
}

//
//.. AddIndex
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::AddIndex(int i, TIndex new_index)
{
  TMDimIndex<TIndex> I(DIM);
  sub_index[i]->NewEntry() = new_index;
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
    for (size_t i = sub_index[level]->begin_idx();
         i < sub_index[level]->end_idx();
         i = sub_index[level]->next_idx(i)) {
      TMDimIndex<TIndex> next_I = I;
      next_I[level] = sub_index[level]->At(i);
      CreateTotalIndex(level + 1, next_I);
    };
  } else {
    try {
      total_index->Find(I);
    } 
    catch (NotFound_error) {
      total_index->NewEntry(I) = NULL;
    };
    //if (total_index->Find(I) == -1) total_index->NewEntry(I) = NULL;
  };
}

//
//.. AllocNewBlocks
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::AllocNewBlocks()
{
  int num_new_blocks = 0;
  TValue *new_root_pointer = new TValue[NumBlocks() * MaxNumEntries()];

  // find out hom many new blocks are needed
  FORALL(i, total_index->) {
    if (total_index->ViAt(i) == NULL) num_new_blocks++;
  };
                                       
  // copy existing data
  //size_t h;
  for (size_t block = 0; block < NumBlocks() - num_new_blocks; block++) {
    for (size_t i = 0; i < MaxNumEntries(); i++) {
      //h = MAP::TotalIndex(this, block, i);
      new_root_pointer[MAP::TotalIndex(MaxNumEntries(), NumBlocks(), block, i)] =
        root_pointer[MAP::TotalIndex
                    (MaxNumEntries(), NumBlocks() - num_new_blocks, block, i)];
    };
  };

  // set recently allocated data to the default value
  for (size_t block = NumBlocks() - num_new_blocks; block < NumBlocks(); block++) {
    for (size_t i = 0; i < MaxNumEntries(); i++) {
      new_root_pointer[MAP::TotalIndex(this, block, i)] = default_value;
    };
  };

  // take the new pointer
  delete [] root_pointer;
  root_pointer = new_root_pointer;

  // set the entries in total_index to the correct positions in the new field
  FORALL(block, total_index->) {
    total_index->ViAt(block) = root_pointer + MAP::TotalIndex(this, block, 0);
  };

  // update all connected variables
  FORALL(i, mapped_vars->) {
    mapped_vars->At(i)->Update();
  };
}

//
//.. Extend
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::Extend(size_t delta)
{
  TValue* new_root_pointer = new TValue[NumBlocks() * (MaxNumEntries() + delta)];

  // find out how many entries need to be copied
  size_t num_copy_entries = MaxNumEntries();
  if (delta < 0) num_copy_entries -= delta;

  // copy existing data to the new field
  for(size_t block = total_index->begin_idx();
      block < total_index->end_idx();
      block = total_index->next_idx(block)) {
    for (size_t i = 0; i < num_copy_entries; i++) {
      new_root_pointer[MAP::TotalIndex(MaxNumEntries() + delta, NumBlocks(), block, i)] =
        root_pointer[MAP::TotalIndex(MaxNumEntries(), NumBlocks(), block, i)];
    };
  };

  // take the new pointer
  delete [] root_pointer;
  root_pointer = new_root_pointer;

  // set new data to the default value
  for(size_t block = total_index->begin_idx();
      block < total_index->end_idx();
      block = total_index->next_idx(block)) {
    for (size_t i = num_copy_entries; i < MaxNumEntries() + delta; i++) {
      root_pointer[MAP::TotalIndex(MaxNumEntries() + delta, NumBlocks(), block, i)] 
        = default_value;
    };
  };

  // set the entries in total_index to the correct positions in the new field
  for(size_t block = total_index->begin_idx();
      block < total_index->end_idx();
      block = total_index->next_idx(block)) {
    total_index->ViAt(block) = 
      root_pointer + MAP::TotalIndex(MaxNumEntries() + delta, NumBlocks(), block, 0);
  };

  // update all connected variables
  for(size_t i = mapped_vars->begin_idx();
      i < mapped_vars->end_idx();
      i = mapped_vars->next_idx(i)) {
    mapped_vars->At(i)->Update();
  };
}

//
//.. CopyEntry
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::CopyEntry(size_t src, size_t dst)
{
  List::CopyEntry(src, dst);
  FORALL(i, total_index->) {
    root_pointer[i * MaxNumEntries() + dst] = root_pointer[i * MaxNumEntries() + src];
  };
}

//
//.. InitEntry
//
template<class TValue, class TIndex, int DIM, class MAP>
void TMDimList<TValue, TIndex, DIM, MAP>::InitEntry(size_t i)
{
  FORALL(j, total_index->) total_index->ViAt(j)[MAP::RelativeIndex(this, i)] 
    = default_value;
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
    FORALL(i, sub_index[level]->) {
      sub_ind = GetSubIndices(I + sub_index[level]->At(i));
      FORALL(j, sub_ind->) ind->NewEntry() = sub_ind->At(j);
      delete sub_ind;
    };
  };
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
    return total_index->At(I);
  };
}

//
//.. IndexExists
//
template<class TValue, class TIndex, int DIM, class MAP>
bool TMDimList<TValue, TIndex, DIM, MAP>::IndexExists(size_t level, TIndex an_index)
{
  bool found = false;
  FORALL(i, sub_index[level]->) {
    if (sub_index[level]->At(i) == an_index) found = true;
  };
  return found;
}

} // namespace

#endif









