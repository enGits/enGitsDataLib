// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2024 enGits GmbH                                         +
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

#ifndef SHORTDELTAVECTOR_H
#define SHORTDELTAVECTOR_H

#include "edl/edl.h"
#include <bits/stdint-uintn.h>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <utility> // for std::swap
#include <iostream>

/**
  * @brief ShortDeltaVector is a vector like container with very low memory overhead.
  * It is designed to store integer like data (e.g. pointers). Only the difference to
  * a reference value is stored. The reference value will be automatically adapted
  * such that a smaller data type can be used than the target data type.
  * The maximal length of the vector is 65536 (16 bit)
  */
template <class TValue, class TIndex=uint8_t>
class ShortDeltaVector
{
public: // data types

  typedef TValue value_type;
  typedef TIndex index_type;


private: // attributes

  static const index_type m_StartSize = 8;

  index_type m_VectorSize    = 0;
  index_type m_AllocatedSize = 0;
  uint8_t    m_DeltaSize     = 1;
  uint8_t*   m_Data          = nullptr;
  uintptr_t  m_Reference     = 0;


public: // attributes

#ifdef PLATFORM_BIG_ENDIAN
  static const bool m_IsBigEndian = true;
#else
  static const bool m_IsBigEndian = false;
#endif


private: // methods

  size_t entryIndex(index_type i) const
  {
    return i*m_DeltaSize;
  }

  void grow(index_type new_items)
  {
    index_type  old_size  = m_VectorSize;
    uint64_t    new_size  = old_size + new_items;
    uint8_t*    new_data  = new uint8_t[new_size*m_DeltaSize];
    size_t      loop_size = m_AllocatedSize*m_DeltaSize;
    //
    if (!sizeValid(new_size)) {
      throw std::runtime_error("ShortDeltaVector: new size exceeds maximum size");
    }
    //
    for (size_t i = 0; i < loop_size; ++i) {
      new_data[i] = m_Data[i];
    }
    delete[] m_Data;    
    m_Data = new_data;
    m_AllocatedSize = new_size;
  }

  void changeDeltaSize(uint8_t new_delta_size)
  {
    uint8_t* new_data = new uint8_t[m_AllocatedSize*new_delta_size];
    for (size_t i = 0; i < m_VectorSize; ++i) {
      auto delta = get(&m_Data[entryIndex(i)], m_DeltaSize);
      rawSet(delta, &new_data[i * new_delta_size], new_delta_size);
    }
    delete[] m_Data;
    m_Data = new_data;
    m_DeltaSize = new_delta_size;
  }

  void updateReference(value_type value)
  {
    uintptr_t value_ptr = reinterpret_cast<uintptr_t>(value);
    if (m_VectorSize == 0) {
      m_Reference = value_ptr;
    } else {
      uint64_t max_entry     = 0;
      uint64_t delta         = 0;
      bool     update_deltas = false;
      //
      if (value_ptr < m_Reference) {
        update_deltas = true;
        delta         = m_Reference - value_ptr;
        m_Reference   = value_ptr;
        //
        for (index_type i = 0; i < m_VectorSize; ++i) {
          auto idx = entryIndex(i);
          max_entry = std::max(get(&m_Data[idx], m_DeltaSize) + delta, max_entry);
        }
      } else {
        max_entry = value_ptr - m_Reference;
      }
      //
      uint64_t max_delta = maxDelta(m_DeltaSize);
      if (max_entry > max_delta) {
        uint8_t new_delta_size = m_DeltaSize;
        for (uint8_t i = m_DeltaSize + 1; i <= 8; ++i) {
          max_delta = maxDelta(i);
          // if (i > 5) {
          //   std::cout << "**** " << int(i) << " : " << value << ", " << m_Reference << ", " << max_entry << ", " << value - m_Reference << std::endl;
          //   int dummy;
          //   std::cin >> dummy;
          //   std::cout << "**** " << int(i) << " : " << value << ", " << m_Reference << ", " << max_entry << ", " << value - m_Reference << std::endl;
          //   std::cout << std::endl;
          // }
          if (max_entry <= max_delta) {
            new_delta_size = i;
            break;
          }
        }
        changeDeltaSize(new_delta_size);
      }
      //
      if (update_deltas) {
        for (index_type i = 0; i < m_VectorSize; ++i) {
          auto idx = entryIndex(i);
          auto new_value = get(&m_Data[idx], m_DeltaSize) + delta;
          rawSet(new_value, &m_Data[idx], m_DeltaSize);
        }
      }
    }
  }


public: // methods

  static bool sizeValid(size_t desired_size)
  {
    return desired_size <= std::numeric_limits<index_type>::max();
  }

  static uint64_t get(const uint8_t* data, const uint8_t length)
  {
    uint64_t value(0);
#ifdef PLATFORM_BIG_ENDIAN
    for (size_t i = 0; i < length; ++i) {
      reinterpret_cast<uint8_t*>(&value)[i] = data[sizeof(value_type) - i - 1];
    }
#else
    for (size_t i = 0; i < length; ++i) {
      reinterpret_cast<uint8_t*>(&value)[i] = data[i];
    }
#endif
    return value;
  }

  static void rawSet(uint64_t value, uint8_t *data, const uint8_t length) 
  {
#ifdef PLATFORM_BIG_ENDIAN
    for (size_t i = 0; i < length; ++i) {
      data[sizeof(uint64_t) - i - 1] = reinterpret_cast<uint8_t*>(&value)[i];
    }
#else
    for (size_t i = 0; i < length; ++i) {
      data[i] = reinterpret_cast<uint8_t*>(&value)[i];
    }
#endif
  }

  static uint64_t maxDelta(uint8_t delta_size)
  {
    uint64_t max_value = std::numeric_limits<uint8_t>::max() + 1;
    for (uint8_t i = 1; i < delta_size; ++i) {
      max_value *= std::numeric_limits<uint8_t>::max() + 1;
    }
    return max_value - 1;
  }

  ShortDeltaVector(size_t initial_size = 0)
  {
    index_type alloc_size = initial_size;    
    if (alloc_size == 0) {
      alloc_size = m_StartSize;
    }
    if (!sizeValid(alloc_size)) {
      throw std::runtime_error("ShortDeltaVector: initial size exceeds maximum size");
    }
    if (!sizeValid(initial_size)) {
      throw std::runtime_error("ShortDeltaVector: initial size exceeds maximum size");
    }
    m_Data = new std::uint8_t[alloc_size*m_DeltaSize];
    m_AllocatedSize = alloc_size;
  }

  ~ShortDeltaVector()
  {
    delete[] m_Data;
  }

  uint16_t deltaSize() const
  {
    return m_DeltaSize;
  }

  bool checkDeltaSize() const
  {
    uint64_t min_delta = std::numeric_limits<uint64_t>::max();
    uint64_t max_delta = 0;
    for (index_type i = 0; i < m_VectorSize; ++i) {
      auto idx = entryIndex(i);
      auto delta = get(&m_Data[idx], m_DeltaSize);
      max_delta = std::max(delta, max_delta);
      min_delta = std::min(delta, min_delta);
    }
    if (min_delta != 0) {
      return false;
    }
    uint8_t expected_delta_size = 0;
    for (uint8_t i = 1; i <= 8; ++i) {
      if (max_delta <= maxDelta(i)) {
        expected_delta_size = i;
        break;
      }
    }
    if (m_DeltaSize != expected_delta_size) {
      return false;
    }
    return true;
  }

  void resize(size_t new_size)
  {
    if (!sizeValid(new_size)) {
      throw std::runtime_error("ShortDeltaVector: new size exceeds maximum size");
    }
    index_type loop_size = std::min(m_VectorSize, index_type(new_size));
    uint8_t*   new_data  = new uint8_t[new_size*m_DeltaSize];
    for (size_t i = 0; i < loop_size*m_DeltaSize; ++i) {
      new_data[i] = m_Data[i];
    }
    delete[] m_Data;
    m_Data = new_data;
    m_AllocatedSize = new_size;
    m_VectorSize    = new_size;
  }

  void reserve(size_t new_size)
  {
    if (!sizeValid(new_size)) {
      throw std::runtime_error("ShortDeltaVector: new size exceeds maximum size");
    }
    index_type loop_size = std::min(m_VectorSize, index_type(new_size));
    uint8_t*   new_data  = new uint8_t[new_size*m_DeltaSize];
    for (size_t i = 0; i < loop_size*m_DeltaSize; ++i) {
      new_data[i] = m_Data[i];
    }
    delete[] m_Data;
    m_Data = new_data;
    m_AllocatedSize = new_size;
    m_VectorSize    = loop_size;
  }

  void clear()
  {
    m_VectorSize = 0;
  }

  size_t size() const
  {
    return m_VectorSize;
  }

  void push_back(const value_type& value)
  {
    if (m_VectorSize == m_AllocatedSize) {
      grow(m_AllocatedSize);
    }
    updateReference(value);
    auto idx = entryIndex(m_VectorSize);
    uintptr_t value_ptr = reinterpret_cast<uintptr_t>(value);
    uint64_t  delta = value_ptr - m_Reference;
    rawSet(delta, &m_Data[idx], m_DeltaSize);
    ++m_VectorSize;
  }

  value_type operator[](size_t i)
  {
    if (m_VectorSize != m_AllocatedSize) {
      shrink_to_fit();
    }
    auto idx = entryIndex(i);
    uint64_t delta = get(&m_Data[idx], m_DeltaSize);
    return reinterpret_cast<value_type>(m_Reference + delta);
  }

  bool operator==(const ShortDeltaVector& other) const
  {
    if (m_VectorSize != other.m_VectorSize) {
      return false;
    }
    for (index_type i = 0; i < m_VectorSize; ++i) {
      auto idx = entryIndex(i);
      auto value1 = get(&m_Data[idx], m_DeltaSize) + m_Reference;
      auto value2 = get(&other.m_Data[idx], m_DeltaSize) + other.m_Reference;
      if (value1 != value2) {
        return false;
      }
    }
    return true;
  }

  bool operator!=(const ShortDeltaVector& other) const
  {
    return !(*this == other);
  }

  void set(size_t i, const value_type& value)
  {
    if (i >= m_VectorSize) {
      throw std::runtime_error("ShortDeltaVector: index out of range");
    }
    updateReference(value);
    auto idx = entryIndex(i);
    rawSet(value - m_Reference, &m_Data[idx], m_DeltaSize);
  }

  value_type operator[](size_t i) const
  {
    auto idx = entryIndex(i);
    uint64_t delta = get(&m_Data[idx], m_DeltaSize);
    return reinterpret_cast<value_type>(m_Reference + delta);
  }

  void shrink_to_fit() 
  {
    index_type new_size = m_VectorSize;
    uint8_t* new_data = new uint8_t[new_size*m_DeltaSize];
    for (size_t i = 0; i < new_size*m_DeltaSize; ++i) {
      new_data[i] = m_Data[i];
    }
    delete[] m_Data;
    m_Data = new_data;
    m_AllocatedSize = new_size;
  }
  
  void swap(index_type i1, index_type i2)
  {
    if (i1 < 0 || i1 >= m_VectorSize || i2 < 0 || i2 >= m_VectorSize) {
      throw std::runtime_error("ShortDeltaVector: index out of range");
    }
    auto idx1 = entryIndex(i1);
    auto idx2 = entryIndex(i2);
    for (uint8_t i = 0; i < m_DeltaSize; ++i) {
      std::swap(m_Data[idx1 + i], m_Data[idx2 + i]);
    }
  }

  size_t memoryUsage() const
  {
    return m_AllocatedSize*m_DeltaSize + sizeof(*this);
  }

  void debugPrint() const
  {
    using namespace std;
    cout << endl;
    cout << "size           = " << int(m_VectorSize) << endl;
    cout << "allocated size = " << int(m_AllocatedSize) << endl;
    cout << "reference      = " << m_Reference << endl;
    cout << "delta size     = " << int(m_DeltaSize) << endl;
    for (index_type i = 0; i < m_VectorSize; ++i) {
      auto idx = entryIndex(i);
      auto delta = get(&m_Data[idx], m_DeltaSize);
      cout << "  " << i << " : " << delta << "/" << delta + m_Reference << endl;
    }
  }


public: // iterators

  struct const_iterator
  {
    size_t i;
    const ShortDeltaVector<value_type,index_type>* dv;
    void operator++() { ++i; }
    bool operator!=(const const_iterator& other) const 
    { 
      return (i != other.i) || (dv != other.dv); 
    }
    value_type operator*() const 
    { 
      return dv->operator[](i); 
    }
  };

  const_iterator begin() const
  {
    return {0, this};
  }

  const_iterator end() const
  {
    return {m_VectorSize, this};
  }

  const_iterator cbegin() const
  {
    return {0, this};
  }

  const_iterator cend() const
  {
    return {m_VectorSize, this};
  }

};

// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

#include <iostream>
#include <random>
#include "edl/quicktimer.h"

TEST_CASE ("ShortDeltaVector_endianness")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  typedef ShortDeltaVector<uint64_t,uint16_t> deltavec_t;
  size_t N = 100;
  size_t n = 16;
  //
  // create a std::vector with N random values of n bit integers
  //
  vector<uint64_t> data(N);
  uint64_t upper_limit = 2;
  for (size_t i = 1; i < n; ++i) {
    upper_limit *= 2;
  }
  upper_limit -= 1;
  random_device rd;
  mt19937_64 eng(rd());
  uniform_int_distribution<uint64_t> distr(0, upper_limit);
  for (size_t i = 0; i < N; ++i) {
    data[i] = distr(eng);
  }    
  //
  // check that all values can be passed through the static set/get methods of deltavec_t
  //
  for (size_t i = 0; i < N; ++i) {
    uint8_t data_bytes[4];
    deltavec_t::rawSet(data[i], data_bytes, 4);
    uint64_t value = deltavec_t::get(data_bytes, 4);
    CHECK(value == data[i]);
  }
}

TEST_CASE("ShortDeltaVector_memory_overhead") 
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  CHECK(sizeof(ShortDeltaVector<uint64_t,uint8_t>)  == 24);
  CHECK(sizeof(ShortDeltaVector<uint64_t,uint16_t>) == 24);
  CHECK(sizeof(ShortDeltaVector<uint64_t,uint32_t>) == 32);
  CHECK(sizeof(ShortDeltaVector<uint64_t,uint64_t>) == 40);
  //
  CHECK(sizeof(ShortDeltaVector<uint32_t,uint8_t>)  == sizeof(vector<uint64_t>));
}

TEST_CASE("ShortDeltaVector_simple")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  typedef ShortDeltaVector<uint64_t,uint8_t> deltavec_t;
  vector<uint64_t> data = {750, 500, 750, 250, 300, 400, 500, 600, 700, 800, 900};
  //
  // copy the data to a ShortDeltaVector
  //
  deltavec_t dv;
  uint64_t min_value = numeric_limits<uint64_t>::max();
  uint64_t max_value = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    dv.push_back(data[i]);
    min_value = min(min_value, data[i]);
    max_value = max(max_value, data[i]);
  }
  //
  // check if the data is correct
  //
  for (size_t i = 0; i < dv.size(); ++i) {
    CHECK(dv[i] == data[i]);
  }
  //
  // check memory overhead
  //
  CHECK(sizeof(deltavec_t) == 24);
  //
  // compute expected delta size
  //
  uint8_t expected_delta_size = 0;
  uint64_t max_delta = max_value - min_value;
  for (uint8_t i = 1; i < 8; ++i) {
    if (max_delta <= (1 << (8*i)) - 1) {
      expected_delta_size = i;
      break;
    }
  }
  //
  // check if the memory usage is correct
  //    
  size_t expected_memory = expected_delta_size*data.size() + sizeof(deltavec_t);
  CHECK(expected_delta_size == dv.deltaSize());
  CHECK(dv.memoryUsage() == expected_memory);
}

TEST_CASE("ShortDeltaVector_simple_large")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  typedef ShortDeltaVector<uint64_t,uint16_t> deltavec_t;
  size_t N = 10000;
  //
  for (size_t n = 4; n <= 64; n += 4) {
    //
    // create a std::vector with N random values of n bit integers stored as 64 bit integers
    //
    vector<uint64_t> data(N);
    uint64_t upper_limit = 2;
    for (size_t i = 1; i < n; ++i) {
      upper_limit *= 2;
    }
    upper_limit -= 1;
    random_device rd;
    mt19937_64 eng(rd());
    uniform_int_distribution<uint64_t> distr(0, upper_limit);
    for (size_t i = 0; i < N; ++i) {
      data[i] = distr(eng);
      CHECK(data[i] <= upper_limit);
    }    
    //
    // copy the data to a ShortDeltaVector and compute the maximal and minimal values
    //
    deltavec_t dv;
    uint64_t min_value = numeric_limits<uint64_t>::max();
    uint64_t max_value = 0;
    for (size_t i = 0; i < N; ++i) {
      dv.push_back(data[i]);
      min_value = min(min_value, data[i]);
      max_value = max(max_value, data[i]);
    }
    //
    // check if the data is correct
    //
    for (size_t i = 0; i < N; ++i) {
      CHECK(dv[i] == data[i]);
    }
    //
    // compute expected delta size
    //
    uint8_t expected_delta_size = 0;
    uint64_t max_delta = max_value - min_value;
    for (uint8_t i = 1; i <= 8; ++i) {
      if (max_delta <= dv.maxDelta(i)) {
        expected_delta_size = i;
        break;
      }
    }
    //
    // check if the memory usage is correct
    //    
    size_t expected_memory = expected_delta_size*N + sizeof(deltavec_t);
    CHECK(expected_delta_size == dv.deltaSize());
    CHECK(dv.memoryUsage() == expected_memory);
  }
}

TEST_CASE("ShortDeltaVector_access_performance")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  typedef double real_t;
  typedef ShortDeltaVector<real_t*,uint16_t> deltavec_t;
#ifdef EDL_DEBUG
  size_t N = 10000;
  size_t number_of_loops = 1;
#else 
  size_t N = 10000;
  size_t number_of_loops = 1000;
#endif
  //
  // create a std::vector with N pointers to random real_t values between 0 and 1
  //
  vector<real_t*> data(N);
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dist(0.0, 1.0);
  real_t* max_entry(0);
  real_t* min_entry = max_entry - 1;
  for (size_t i = 0; i < N; ++i) {
    data[i] = new real_t(dist(gen));
    max_entry = max(max_entry, (data[i]));
    min_entry = min(min_entry, (data[i]));
  }
  cout << "max entry : " << max_entry << endl;
  cout << "min entry : " << min_entry << endl;
  cout << "max delta : " << max_entry - min_entry << endl;
  //
  // copy the pointers to a ShortDeltaVector
  //
  deltavec_t dv;
  real_t* reference = nullptr;
  for (size_t i = 0; i < N; ++i) {
    if (i == 0) {
      reference = data[i];
    } else {
      reference = min(reference, data[i]);
    }
    dv.push_back(data[i]);
  }
  CHECK(max_entry - min_entry <= dv.maxDelta(dv.deltaSize()));
  //
  // perform several loops of square and square root operations and measure the time
  //   - make a copy of the original data
  //   - first on the pointers in the std::vector
  //   - second on the pointers in the ShortDeltaVector
  //   - always compare the results to the original values
  //
  vector<real> reference_data(N);
  for (size_t i = 0; i < N; ++i) {
    reference_data[i] = *data[i];
  }
  //
  // compute on vector<real*>
  //
  QuickTimer timer;
  timer.start();
  for (size_t i_loop = 0; i_loop < number_of_loops; ++i_loop) {
    for (size_t i = 0; i < N; ++i) {
      *data[i] = (*data[i])*(*data[i]);
    }
    for (size_t i = 0; i < N; ++i) {
      *data[i] = sqrt(*data[i]);
    }    
  }
  timer.stop();
  for (size_t i = 0; i < N; ++i) {
    CHECK(abs(reference_data[i] - *data[i]) < 1e-20);
  }
  auto msecs1 = timer.milliseconds();
  timer.reset();
  //
  // compute on ShortDeltaVector<real*>
  //
  timer.start();
  for (size_t i_loop = 0; i_loop < number_of_loops; ++i_loop) {
    for (size_t i = 0; i < N; ++i) {
      *dv[i] = (*dv[i])*(*dv[i]);
    }
    for (size_t i = 0; i < N; ++i) {
      *dv[i] = sqrt(*dv[i]);
    }    
  }
  timer.stop();   
  for (size_t i = 0; i < N; ++i) {
    CHECK(abs(reference_data[i] - *dv[i]) < 1e-20);
  }
  auto msecs2 = timer.milliseconds();
  //
  cout << "delta size                : " << dv.deltaSize() << endl;
  cout << "number of accesses        : " << N*number_of_loops << endl;
  cout << "time for std::vector      : " << msecs1 << " ms" << endl;
  cout << "time for ShortDeltaVector : " << msecs2 << " ms" << endl;
  //
  // delete the data
  //
  for (size_t i = 0; i < N; ++i) {
    delete data[i];
  }
}

TEST_CASE("ShortDeltaVector_iterators")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  typedef ShortDeltaVector<uint64_t,uint16_t> deltavec_t;
  size_t N = 10000;
  size_t n = 32;
  //
  // create a std::vector with N random values of n bit integers stored as 64 bit integers
  //
  vector<uint64_t> data(N);
  uint64_t upper_limit = 2;
  for (size_t i = 1; i < n; ++i) {
    upper_limit *= 2;
  }
  upper_limit -= 1;
  random_device rd;
  mt19937_64 eng(rd());
  uniform_int_distribution<uint64_t> distr(0, upper_limit);
  for (size_t i = 0; i < N; ++i) {
    data[i] = distr(eng);
    CHECK(data[i] <= upper_limit);
  }    
  //
  // copy the data to a ShortDeltaVector and compute the maximal and minimal values
  //
  deltavec_t dv;
  uint64_t min_value = numeric_limits<uint64_t>::max();
  uint64_t max_value = 0;
  for (size_t i = 0; i < N; ++i) {
    dv.push_back(data[i]);
    min_value = min(min_value, data[i]);
    max_value = max(max_value, data[i]);
  }
  //
  // check if the data is correct (range iteration)
  //
  {
    size_t i = 0;
    for (auto v : dv) {
      CHECK(v == data[i]);
      ++i;
    }
  }
  //
  // check if the data is correct (iterator)
  //
  {
    size_t i = 0;
    for (auto it = dv.begin(); it != dv.end(); ++it) {
      CHECK(*it == data[i]);
      ++i;
    }
  }
}

TEST_CASE("ShortDeltaVector_set")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  typedef ShortDeltaVector<uint64_t,uint16_t> deltavec_t;
  //
  vector<uint64_t> data = {750, 500, 750, 250, 300, 400, 500, 600, 700, 800};
  //
  // copy the data to a ShortDeltaVector and compute the maximal and minimal values
  //
  deltavec_t dv;
  for (size_t i = 0; i < data.size(); ++i) {
    dv.push_back(data[i]);
  }
  //
  // loop over data and increment each value by 2**16
  //
  for (size_t i = 0; i < data.size(); ++i) {
    data[i] += uint64_t(1) << 16;
    dv.set(i, data[i]);
  }
  //
  // check if the data is correct
  //
  for (size_t i = 0; i < data.size(); ++i) {
    CHECK(dv[i] == data[i]);
  }
  //
  // loop over data and increment each value by 2**32
  //
  for (size_t i = 0; i < data.size(); ++i) {
    data[i] += uint64_t(1) << 32;
    dv.set(i, data[i]);
  }
  //
  // check if the data is correct
  //
  for (size_t i = 0; i < data.size(); ++i) {
    CHECK(dv[i] == data[i]);
  }
}

TEST_CASE("ShortDeltaVector_clear_and_comparison")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  typedef ShortDeltaVector<uint64_t,uint16_t> deltavec_t;
  //
  vector<uint64_t> data1 = {750, 500, 750, 250, 300, 400, 500, 600, 700, 800};
  vector<uint64_t> data2 = {800, 550, 800, 300, 350, 450, 550, 650, 750, 850};
  //
  // copy the data to a ShortDeltaVector and compute the maximal and minimal values
  //
  deltavec_t dv1, dv2;
  for (size_t i = 0; i < data1.size(); ++i) {
    dv1.push_back(data1[i]);
    dv2.push_back(data2[i]);
  }
  CHECK(dv1 != dv2);
  for (size_t i = 0; i < data1.size(); ++i) {
    auto v = dv1[i];
    v += 50;
    dv1.set(i, v);
  }
  CHECK(dv1 == dv2);
  dv1.clear();
  CHECK(dv1 != dv2);
  for (size_t i = 0; i < data1.size(); ++i) {
    dv1.push_back(data1[i]);
  }
  CHECK(dv1 != dv2);
  for (size_t i = 0; i < data1.size(); ++i) {
    auto v = dv1[i];
    v += 50;
    dv1.set(i, v);
  }
  CHECK(dv1 == dv2);
}

/*
TEST_CASE("ShortDeltaVector_find")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  typedef ShortDeltaVector<uint64_t,uint16_t> deltavec_t;
  //
  vector<uint64_t> data = {750, 500, 750, 250, 300, 400, 500, 600, 700, 800};
  //
  // copy the data to a ShortDeltaVector and compute the maximal and minimal values
  //
  deltavec_t dv;
  for (size_t i = 0; i < data.size(); ++i) {
    dv.push_back(data[i]);
  }
  //
  // find the value 500
  //
  auto it = find(dv.begin(), dv.end(), 500);
  CHECK(it != dv.end());
  CHECK(*it == 500);  
}
*/

TEST_CASE("ShortDeltaVector_swap")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  typedef ShortDeltaVector<uint64_t,uint16_t> deltavec_t;
  //
  vector<uint64_t> data = {750, 500, 750, 250, 300, 400, 500, 600, 700, 800};
  //
  // copy the data to a ShortDeltaVector
  //
  deltavec_t dv;
  for (size_t i = 0; i < data.size(); ++i) {
    dv.push_back(data[i]);
  }
  //
  // check swap
  //
  CHECK(dv[3] == 250);
  CHECK(dv[3] != dv[4]);
  dv.swap(3, 4);
  CHECK(dv[3] == 300);
  CHECK(dv[4] == 250);
}

#endif // SHORTDELTAVECTOR_H