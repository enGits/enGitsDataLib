// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2016 enGits GmbH                                         +
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

#ifndef SHORTVECTOR_H
#define SHORTVECTOR_H

#include "edl/edl.h"
#include "edl/mathvector.h"
#include <bits/stdint-uintn.h>
#include <cstddef>
#include <cstdint>

/**
  * @brief ShortVector is a vector like container with very low memory overhead.
  *
  * entry 0 : current size of the vector
  * entry 1 : current size of the allocated array
  * entry 2 : first element
  * entry 3 : second element
  * ...
  */
namespace EDL_NAMESPACE
{

template <class TValue, class TEntry=std::uint8_t>
class ShortVector
{

public: // data types

  typedef TValue value_type;
  typedef TEntry entry_type;

  typedef value_type* iterator;
  typedef const value_type* const_iterator;


private: // attributes
  
  entry_type* m_Data;
  
  static const std::uint8_t m_ItemSize   = sizeof(TValue)/sizeof(TEntry);
  static const std::uint8_t m_EntrySize  = sizeof(TEntry);
  static const std::uint8_t m_StartSize  = 8;
  static const std::uint8_t m_HeaderSize = 2;


private: // methods

  entry_type allocatedSize() const
  {
    return m_Data[1];
  }

  entry_type& allocatedSize()
  {
    return m_Data[1];
  }

  entry_type arraySize() const
  {
    return allocatedSize()*m_ItemSize + m_HeaderSize;
  }

  entry_type& vectorSize()
  {
    return m_Data[0];
  }

  entry_type vectorSize() const
  {
    return m_Data[0];
  }

  entry_type entryIndex(entry_type i) const
  {
    return i*m_ItemSize  + m_HeaderSize;
  }

  void grow(entry_type new_items)
  {
    entry_type  old_size  = vectorSize();
    entry_type  new_size  = old_size + new_items;
    entry_type* new_data  = new entry_type[new_size*m_ItemSize + m_HeaderSize];
    entry_type  loop_size = allocatedSize()*m_ItemSize + m_HeaderSize;
    //
    if (!sizeValid(new_size)) {
      throw std::runtime_error("ShortVector: new size exceeds maximum size");
    }
    //
    for (entry_type i = 0; i < loop_size; ++i) {
      new_data[i] = m_Data[i];
    }
    delete[] m_Data;    
    m_Data = new_data;
    allocatedSize() = new_size;
  }


public: // methods

  ShortVector(size_t initial_size = 0)
  {
    entry_type alloc_size = initial_size;    
    if (alloc_size == 0) {
      alloc_size = m_StartSize;
    }
    if (!sizeValid(alloc_size)) {
      throw std::runtime_error("ShortVector: initial size exceeds maximum size");
    }
    if (!sizeValid(initial_size)) {
      throw std::runtime_error("ShortVector: initial size exceeds maximum size");
    }
    m_Data = new entry_type[alloc_size*m_ItemSize  + m_HeaderSize];
    m_Data[0] = initial_size;
    m_Data[1] = alloc_size;
  }

  ~ShortVector()
  {
    delete[] m_Data;
  }

  void resize(size_t new_size)
  {
    if (!sizeValid(new_size)) {
      throw std::runtime_error("ShortVector: new size exceeds maximum size");
    }
    entry_type  loop_size = std::min(vectorSize(), entry_type(new_size));
    entry_type* new_data = new entry_type[new_size*m_ItemSize + m_HeaderSize];
    for (entry_type i = 0; i < loop_size*m_ItemSize + m_HeaderSize; ++i) {
      new_data[i] = m_Data[i];
    }
    delete[] m_Data;
    m_Data = new_data;
    allocatedSize() = new_size;
    vectorSize()    = new_size;
  }

  size_t size() const
  {
    return m_Data[0];
  }

  void push_back(const value_type& value)
  {
    if (vectorSize() == allocatedSize()) {
      grow(allocatedSize());
    }
    auto idx = entryIndex(vectorSize());
    value_type* ptr = reinterpret_cast<value_type*>(&m_Data[idx]);
    ++vectorSize();
    *ptr = value;
  }

  value_type& operator[](entry_type i)
  {
    if (vectorSize() != allocatedSize()) {
      compress();
    }
    value_type* ptr = reinterpret_cast<value_type*>(&m_Data[entryIndex(i)]);
    return *ptr;
  }

  value_type operator[](entry_type i) const
  {
    if (vectorSize() != allocatedSize()) {
      compress();
    }
    value_type* ptr = reinterpret_cast<value_type*>(&m_Data[entryIndex(i)]);
    return *ptr;
  }

  void compress()
  {
    entry_type new_size = vectorSize();
    entry_type* new_data = new entry_type[new_size*m_ItemSize + m_HeaderSize];
    for (entry_type i = 0; i < new_size*m_ItemSize + m_HeaderSize; ++i) {
      new_data[i] = m_Data[i];
    }
    delete[] m_Data;
    m_Data = new_data;
    allocatedSize() = new_size;
  }

  size_t memoryUsage() const
  {
    return arraySize()*m_EntrySize + sizeof(entry_type*);
  }

  static bool sizeValid(size_t desired_size)
  {
    size_t max_size = (std::numeric_limits<entry_type>::max() - 2)/m_ItemSize;
    return desired_size <= max_size;
  }

  iterator begin()
  {
    return reinterpret_cast<value_type*>(&m_Data[m_HeaderSize]);
  }

  const_iterator begin() const
  {
    return reinterpret_cast<const value_type*>(&m_Data[m_HeaderSize]);
  }

  const_iterator cbegin() const
  {
    return reinterpret_cast<const value_type*>(&m_Data[m_HeaderSize]);
  }

  iterator end()
  {
    return reinterpret_cast<value_type*>(&m_Data[entryIndex(vectorSize())]);
  }

  const_iterator end() const
  {
    return reinterpret_cast<const value_type*>(&m_Data[entryIndex(vectorSize())]);
  }

  const_iterator cend() const
  {
    return reinterpret_cast<const value_type*>(&m_Data[entryIndex(vectorSize())]);
  }

};


} // namespace EDL_NAMESPACE

// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

#include <iostream>
#include <random>

TEST_CASE("ShortVector_simple")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  CHECK(sizeof(ShortVector<int>) == 8);
  //
  typedef ShortVector<int> shortvec_t;
  shortvec_t sv;
  CHECK(sv.size() == 0);
  int N = 20;
  for (int i = 0; i < N; ++i) {
    sv.push_back(i);
  }
  CHECK(sv.size() == N);
  for (int i = 0; i < N; ++i) {
    CHECK(sv[i] == i);
  }
  auto expected_memory = sizeof(void*) + N*sizeof(shortvec_t::value_type) + 2*sizeof(shortvec_t::entry_type);
  CHECK(sv.memoryUsage() == expected_memory);
}

TEST_CASE("ShortVector_item_loop")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  // random data
  //
  typedef ShortVector<pair<double,int>, uint16_t> shortvec_t;
  int N = 100;
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dist(-1.0, 1.0);
  vector<shortvec_t::value_type> data(N);
  for (int i = 0; i < N; ++i) {
    data[i] = make_pair(dist(gen), i);
  }
  //
  shortvec_t sv(N);
  for (int i = 0; i < N; ++i) {
    sv[i] = data[i];
  }
  CHECK(sv.size() == N);
  for (int i = 0; i < N; ++i) {
    CHECK(sv[i].first == data[i].first);
    CHECK(sv[i].second == i);
  }
  auto expected_memory = sizeof(void*) + N*sizeof(shortvec_t::value_type) + 2*sizeof(shortvec_t::entry_type);
  CHECK(sv.memoryUsage() == expected_memory);
  //
  // item loop
  //
  for (auto& item : sv) {
    int  i = item.second;
    auto x = item.first;
    CHECK(data[i].first == x);
  }
}

TEST_CASE("ShortVector_large_random")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  typedef double real;
  typedef MathVector<StaticVector<real,3>> vec_t;
  int N = 50000;
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dist(-1.0, 1.0);
  //
  // random test data
  //
  vector<vec_t> data;
  for (int i = 0; i < N; ++i) {
    vec_t v;
    for (int j = 0; j < 3; ++j) {
      v[j] = dist(gen);
    }
    data.push_back(v);
  }
  //
  // populate not so short ShortVector with random data
  //
  typedef ShortVector<vec_t, uint32_t> shortvec_t;
  shortvec_t sv;
  sv.resize(N);
  for (int i = 0; i < N; ++i) {
    sv[i] = data[i];
  }
  CHECK(sv.size() == N);
  for (int i = 0; i < N; ++i) {
    vec3_t x1 = sv[i];
    vec3_t x2 = data[i];
    CHECK((x2 - x1).abs() < 1.0e-10);
  }
  auto expected_memory = sizeof(void*) + N*sizeof(vec_t) + 2*sizeof(shortvec_t::entry_type);
  CHECK(sv.memoryUsage() == expected_memory);
}

#endif // SHORTVECTOR_H