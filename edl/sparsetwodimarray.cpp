// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "sparsetwodimarray.h"

namespace EDL_NAMESPACE
{

SparseTwoDimArray::SparseTwoDimArray(size_t mne, size_t delta) 
  : List(mne, delta)
{
  setup();
  m_MasterArray = this;
}

SparseTwoDimArray::SparseTwoDimArray(List *a_master) 
  : List(a_master)
{
  setup();
  m_MasterArray = this;
}

SparseTwoDimArray::SparseTwoDimArray(SparseTwoDimArray *a_master_array) 
  : List(a_master_array->master())
{
  m_MasterArray = a_master_array;
  m_Count = NULL;
  m_First = NULL;
  m_DummyData = new List(m_MasterArray->m_DummyData);
}

void SparseTwoDimArray::linkSparseTwoDimArray(SparseTwoDimArray *a_master_array, std::string link_name)
{
  link(a_master_array, link_name);
  m_MasterArray = a_master_array;
  m_Count = NULL;
  m_First = NULL;
  m_DummyData = new List(m_MasterArray->m_DummyData);
}

SparseTwoDimArray::~SparseTwoDimArray()
{
  delete m_DummyData;
  delete m_First;
  delete m_Count;
}

void SparseTwoDimArray::setup()
{
  m_Count = new TList<size_t>(this, 0);
  m_First = new TList<size_t>(this, 0);
  m_DummyData = new List(maxNumEntries(), maxNumEntries());
}

size_t SparseTwoDimArray::addSubEntry(size_t i)
{
  if (this == m_MasterArray) {
    //size_t j;
    if (!isActive(i)) {
      throw InvalidIndex_error(i);
    } else {
      if ((*m_Count)[i] == 0) {
        //
        // allocate the first entry of this row
        size_t i_data = m_DummyData->alloc(1);
        (*m_First)[i] = i_data;
        (*m_Count)[i]++;
        //
      } else {
        size_t i_data = m_DummyData->tryAlloc(1);
        if (i_data == (*m_First)[i] + (*m_Count)[i]) {
          //
          // block is continuous
          // => add coefficient in a straight forward manner
          i_data = m_DummyData->alloc(1);
          (*m_Count)[i]++;
        } else {
          //
          // block would not be continuous
          // => copy whole block to new location
          i_data = m_DummyData->alloc((*m_Count)[i]+1);
          for (size_t j_data = 0; j_data < (*m_Count)[i]; j_data++) {
            m_DummyData->copy((*m_First)[i] + j_data, i_data + j_data);
          }
          (*m_First)[i] = i_data;
          (*m_Count)[i]++;
        }
      }
      return (*m_Count)[i] - 1;
    }
  } else {
    return m_MasterArray->addSubEntry(i);
  }
}

List* SparseTwoDimArray::dummyData()
{
  if (this == m_MasterArray) {
    return m_DummyData;
  } else {
    return m_MasterArray->m_DummyData;
  }
}

void SparseTwoDimArray::delSubList(size_t i)
{
  if (this == m_MasterArray) {
    m_DummyData->deleteEntries((*m_First)[i], (*m_Count)[i]);
    (*m_First)[i] = 0;
    (*m_Count)[i] = 0;
  } else {
    m_MasterArray->delSubList(i);
  }
}

void SparseTwoDimArray::resize(size_t i, size_t l)
{
  delSubList(i);
  while (l > 0) {
    addSubEntry(i);
    l--;
  }
}

void SparseTwoDimArray::deleteData(size_t i)
{ 
  List::deleteData(i);
  delSubList(i);
}

void SparseTwoDimArray::cleanUp()
{
  if (this == m_MasterArray) {
    List::cleanUp();
    m_DummyData->cleanUp();
  } else {
    m_MasterArray->cleanUp();
  }
}
  
void SparseTwoDimArray::resetData(size_t mne, size_t delta)
{
  if (delta == 0) {
    delta = mne/10;
  }
  m_DummyData->reset(mne,delta);
}

void SparseTwoDimArray::relativeResetData(real mne, real delta)
{
  resetData(size_t(mne*maxNumEntries()),size_t(delta*maxNumEntries()));
}

} // namespace







