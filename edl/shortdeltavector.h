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

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <utility> // for std::swap
#include <iostream>
#include <vector>
#include <cstring>

#include <type_traits>
#if __has_include(<bit>)
  #include <bit>
#endif


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


public: // iterators

  struct const_iterator {
    // --- standard iterator associated types ---
    using iterator_category = std::input_iterator_tag;
    using value_type        = typename ShortDeltaVector::value_type;
    using difference_type   = std::ptrdiff_t;
    using pointer           = void;                           // no operator->; deref returns by value
    using reference         = value_type;                     // by-value is OK for input iterators

    // --- state ---
    const ShortDeltaVector<value_type, index_type>* dv = nullptr;
    size_t i = 0;

    // --- ctors ---
    const_iterator() = default;
    const_iterator(const ShortDeltaVector* d, size_t idx) : dv(d), i(idx) {}

    // --- core ops ---
    reference operator*() const { return (*dv)[i]; }

    const_iterator& operator++() { ++i; return *this; }       // pre-increment
    const_iterator operator++(int) {                          // post-increment
      const_iterator tmp(*this);
      ++(*this);
      return tmp;
    }

    friend bool operator==(const const_iterator& a, const const_iterator& b) {
      return a.dv == b.dv && a.i == b.i;
    }
    friend bool operator!=(const const_iterator& a, const const_iterator& b) {
      return !(a == b);
    }
  };

  const_iterator begin() const
  {
    return const_iterator(this, 0);
  }

  const_iterator end() const
  {
    return const_iterator(this, m_VectorSize);
  }

  const_iterator cbegin() const
  {
    return const_iterator(this, 0);
  }
  
  const_iterator cend() const
  {
    return const_iterator(this, m_VectorSize);
  }


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


  // --- helpers: value <-> integral address ------------------------------------
  static inline uintptr_t toAddr(value_type v) 
  {
    if constexpr (std::is_pointer_v<value_type>) 
    {
      return reinterpret_cast<uintptr_t>(v);
    } else {
      return static_cast<uintptr_t>(v);
    }
  }
  static inline value_type fromAddr(uintptr_t a) 
  {
    if constexpr (std::is_pointer_v<value_type>) {
      return reinterpret_cast<value_type>(a);
    } else {
      return static_cast<value_type>(a);
    }
  }
  // ----------------------------------------------------------------------------

  size_t entryIndex(index_type i) const
  {
    return i*m_DeltaSize;
  }

  void grow(index_type new_items)
  {
    const index_type old_size = m_VectorSize;
    const uint64_t new_size = m_AllocatedSize + new_items; // grow capacity
    if (!sizeValid(new_size)) {
      throw std::runtime_error("ShortDeltaVector: new size exceeds maximum size");
    }
    uint8_t* new_data = new uint8_t[new_size * m_DeltaSize];
    const size_t used_bytes = static_cast<size_t>(old_size) * m_DeltaSize;
    for (size_t i = 0; i < used_bytes; ++i) new_data[i] = m_Data[i];
    delete[] m_Data;
    m_Data = new_data;
    m_AllocatedSize = static_cast<index_type>(new_size);
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
    uintptr_t value_ptr = toAddr(value);
    if (m_VectorSize == 0) {
      m_Reference = value_ptr;
      return;
    }
    uint64_t max_entry = 0;
    uint64_t delta = 0;
    bool update_deltas = false;

    if (value_ptr < m_Reference) {
      update_deltas = true;
      delta = m_Reference - value_ptr;
      m_Reference = value_ptr;
      for (index_type i = 0; i < m_VectorSize; ++i) {
        auto idx = entryIndex(i);
        max_entry = std::max(get(&m_Data[idx], m_DeltaSize) + delta, max_entry);
      }
    } else {
      max_entry = value_ptr - m_Reference;
    }

    uint64_t max_delta = maxDelta(m_DeltaSize);
    if (max_entry > max_delta) {
      uint8_t new_delta_size = m_DeltaSize;
      for (uint8_t i = m_DeltaSize + 1; i <= 8; ++i) {
        if (max_entry <= maxDelta(i)) { new_delta_size = i; break; }
      }
      changeDeltaSize(new_delta_size);
    }

    if (update_deltas) {
      for (index_type i = 0; i < m_VectorSize; ++i) {
        auto idx = entryIndex(i);
        auto new_value = get(&m_Data[idx], m_DeltaSize) + delta;
        rawSet(new_value, &m_Data[idx], m_DeltaSize);
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
    uint64_t value = 0;
#if defined(__cpp_lib_endian) && __cpp_lib_endian >= 201907L
    if (std::endian::native == std::endian::little) {
      for (uint8_t i = 0; i < length; ++i) {
        reinterpret_cast<uint8_t*>(&value)[i] = data[i];
      }
    } else {
      for (uint8_t i = 0; i < length; ++i) {
        reinterpret_cast<uint8_t*>(&value)[length - 1 - i] = data[i];
      }
    }
#else
  #ifdef PLATFORM_BIG_ENDIAN
    for (uint8_t i = 0; i < length; ++i) {
      reinterpret_cast<uint8_t*>(&value)[length - 1 - i] = data[i];
    }
  #else
    for (uint8_t i = 0; i < length; ++i) {
      reinterpret_cast<uint8_t*>(&value)[i] = data[i];
    }
  #endif
#endif
    return value;
  }

  static void rawSet(uint64_t value, uint8_t* data, const uint8_t length)
  {
#if defined(__cpp_lib_endian) && __cpp_lib_endian >= 201907L
    if (std::endian::native == std::endian::little) {
      for (uint8_t i = 0; i < length; ++i) {
        data[i] = reinterpret_cast<uint8_t*>(&value)[i];
      }
    } else {
      for (uint8_t i = 0; i < length; ++i) {
        data[i] = reinterpret_cast<uint8_t*>(&value)[length - 1 - i];
      }
    }
#else
  #ifdef PLATFORM_BIG_ENDIAN
    for (uint8_t i = 0; i < length; ++i) {
      data[i] = reinterpret_cast<uint8_t*>(&value)[length - 1 - i];
    }
  #else
    for (uint8_t i = 0; i < length; ++i) {
      data[i] = reinterpret_cast<uint8_t*>(&value)[i];
    }
  #endif
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

  // Assign from std::vector
  ShortDeltaVector& operator=(const std::vector<value_type>& src) {
    clear();
    reserve(src.size());
    for (const value_type& v : src)
        push_back(v);
    return *this;
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
  
  size_t capacity() const
  {
    return m_AllocatedSize;
  }

  bool empty() const
  {
    return m_VectorSize == 0;
  }

  void push_back(const value_type& value)
  {
    if (m_VectorSize == m_AllocatedSize) {
      grow(m_AllocatedSize);
    }
    updateReference(value);
    auto idx = entryIndex(m_VectorSize);
    uint64_t delta = toAddr(value) - m_Reference;
    rawSet(delta, &m_Data[idx], m_DeltaSize);
    ++m_VectorSize;
  }

  value_type front() const
  {
    return *begin();
  }

  value_type operator[](size_t i)
  {
    auto idx = entryIndex(static_cast<index_type>(i));
    uint64_t delta = get(&m_Data[idx], m_DeltaSize);
    return fromAddr(m_Reference + delta);
  }

  value_type operator[](size_t i) const
  {
    auto idx = entryIndex(static_cast<index_type>(i));
    uint64_t delta = get(&m_Data[idx], m_DeltaSize);
    return fromAddr(m_Reference + delta);
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
    auto idx = entryIndex(static_cast<index_type>(i));
    uint64_t delta = toAddr(value) - m_Reference;
    rawSet(delta, &m_Data[idx], m_DeltaSize);
  }

  // Erase a single element at position 'pos'
  const_iterator erase(const_iterator pos)
  {
    if (pos.dv != this) {
      throw std::runtime_error("ShortDeltaVector::erase: iterator does not belong to this container");
    }
    if (pos.i >= m_VectorSize) {
      throw std::runtime_error("ShortDeltaVector::erase: iterator out of range");
    }

    const index_type idx = static_cast<index_type>(pos.i);

    // Shift the raw bytes to close the gap
    if (idx < m_VectorSize - 1) {
      const size_t start = entryIndex(idx);
      const size_t next  = entryIndex(idx + 1);
      const size_t tail_bytes = (m_VectorSize - idx - 1) * m_DeltaSize;
      std::memmove(&m_Data[start], &m_Data[next], tail_bytes);
    }

    --m_VectorSize;

    // Return iterator to the element that followed the erased one
    return const_iterator(this, idx);
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
    if (i1 >= m_VectorSize || i2 >= m_VectorSize) {
      throw std::runtime_error("ShortDeltaVector: index out of range");
    }
    auto idx1 = entryIndex(i1);
    auto idx2 = entryIndex(i2);
    for (uint8_t i = 0; i < m_DeltaSize; ++i) {
      std::swap(m_Data[idx1 + i], m_Data[idx2 + i]);
    }
  }

  void reverse_inplace() {
    if (m_VectorSize <= 1) return;
    index_type i = 0;
    index_type j = static_cast<index_type>(m_VectorSize - 1);
    while (i < j) {
      swap(i, j);
      ++i;
      --j;
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
    if (max_delta <= (uint64_t{1} << (8*i)) - 1) {
      expected_delta_size = i;
      break;
    }
  }

  //
  // check if the memory usage is correct
  //    
  dv.shrink_to_fit();
  size_t expected_memory = expected_delta_size * data.size() + sizeof(deltavec_t);
  CHECK(dv.deltaSize() == expected_delta_size);
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
    dv.shrink_to_fit();
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

TEST_CASE("ShortDeltaVector_erase")
{
  using namespace EDL_NAMESPACE;
  using namespace std;

  typedef ShortDeltaVector<uint64_t, uint16_t> deltavec_t;

  // Test case 1: Erase a single element from the middle
  {
    vector<uint64_t> data = {10, 20, 30, 40, 50};
    deltavec_t dv;
    dv = data;

    auto it = dv.erase(std::next(dv.begin(), 2)); // Erase the element at index 2 (value 30)

    vector<uint64_t> expected = {10, 20, 40, 50};
    CHECK(dv.size() == expected.size());
    for (size_t i = 0; i < dv.size(); ++i) {
      CHECK(dv[i] == expected[i]);
    }
    CHECK(it == std::next(dv.begin(), 2)); // Iterator should point to the next element (value 40)
    CHECK(*dv.end() == 50);
    CHECK(dv[2] == 40);
  }

  // Test case 2: Erase the first element
  {
    vector<uint64_t> data = {100, 200, 300, 400};
    deltavec_t dv;
    dv = data;

    auto it = dv.erase(dv.begin()); // Erase the first element (value 100)

    vector<uint64_t> expected = {200, 300, 400};
    CHECK(dv.size() == expected.size());
    for (size_t i = 0; i < dv.size(); ++i) {
      CHECK(dv[i] == expected[i]);
    }
    CHECK(it == dv.begin()); // Iterator should point to the new first element (value 200)
  }

  // Test case 3: Erase the last element
  {
    vector<uint64_t> data = {5, 10, 15, 20};
    deltavec_t dv;
    dv = data;

    auto it = dv.erase(std::prev(dv.end())); // Erase the last element (value 20)

    vector<uint64_t> expected = {5, 10, 15};
    CHECK(dv.size() == expected.size());
    for (size_t i = 0; i < dv.size(); ++i) {
      CHECK(dv[i] == expected[i]);
    }
    CHECK(it == dv.end()); // Iterator should point to the end
  }

  // Test case 4: Erase all elements one by one
  {
    vector<uint64_t> data = {1, 2, 3, 4, 5};
    deltavec_t dv;
    dv = data;

    while (!dv.empty()) {
      dv.erase(dv.begin());
    }

    CHECK(dv.size() == 0); // Vector should be empty
  }

  // Test case 5: Erase from an empty vector
  {
    deltavec_t dv;

    CHECK_THROWS_AS(dv.erase(dv.begin()), std::runtime_error); // Should throw an error
  }

  // Test case 6: Erase with an invalid iterator
  {
    vector<uint64_t> data = {10, 20, 30};
    deltavec_t dv;
    dv = data;

    deltavec_t::const_iterator invalid_it;
    CHECK_THROWS_AS(dv.erase(invalid_it), std::runtime_error); // Should throw an error
  }
}

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

TEST_CASE("ShortDeltaVector_reverse_inplace")
{
  using namespace EDL_NAMESPACE;
  using namespace std;

  typedef ShortDeltaVector<uint64_t, uint8_t> deltavec_t;

  // Test case 1: Reverse a vector with multiple elements
  {
    vector<uint64_t> data = {1, 2, 3, 4, 5};
    deltavec_t dv;
    dv = data;

    dv.reverse_inplace();

    vector<uint64_t> expected = {5, 4, 3, 2, 1};
    for (size_t i = 0; i < dv.size(); ++i) {
      CHECK(dv[i] == expected[i]);
    }
  }

  // Test case 2: Reverse a vector with a single element
  {
    vector<uint64_t> data = {42};
    deltavec_t dv;
    dv = data;

    dv.reverse_inplace();

    vector<uint64_t> expected = {42};
    for (size_t i = 0; i < dv.size(); ++i) {
      CHECK(dv[i] == expected[i]);
    }
  }

  // Test case 3: Reverse an empty vector
  {
    deltavec_t dv;

    dv.reverse_inplace();

    CHECK(dv.size() == 0);
  }
}

TEST_CASE("ShortDeltaVector_memory_overhead")
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  cout << " 8 bit index / 64 bit value : " << sizeof(ShortDeltaVector<real*,uint8_t>) << endl;
  cout << "16 bit index / 64 bit value : " << sizeof(ShortDeltaVector<real*,uint16_t>) << endl;
  cout << "32 bit index / 64 bit value : " << sizeof(ShortDeltaVector<real*,uint32_t>) << endl;
  cout << "64 bit index / 64 bit value : " << sizeof(ShortDeltaVector<real*,uint64_t>) << endl;
  cout << " 8 bit index / 32 bit value : " << sizeof(ShortDeltaVector<uint32_t,uint8_t>) << endl;
  cout << "16 bit index / 32 bit value : " << sizeof(ShortDeltaVector<uint32_t,uint16_t>) << endl;
  cout << "32 bit index / 32 bit value : " << sizeof(ShortDeltaVector<uint32_t,uint32_t>) << endl;
  cout << "64 bit index / 32 bit value : " << sizeof(ShortDeltaVector<uint32_t,uint64_t>) << endl;
}

#endif // SHORTDELTAVECTOR_H