// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#include "edl/edl.h"
#include "edl/mathvector.h"
#include <cstddef>
#include <cstdint>
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <algorithm>
#include <limits>
#include <new>
#include <stdexcept>
#include <type_traits>


namespace EDL_NAMESPACE
{

// ---- helpers ---------------------------------------------------------------

template <class T>
constexpr std::size_t alignOf() noexcept { return alignof(T); }

constexpr std::size_t alignUp(std::size_t n, std::size_t a) noexcept {
  return (a == 0) ? n : ((n + (a - 1)) / a) * a;
}

// ---- ShortVector -----------------------------------------------------------

template <class TValue, class TIndex = std::uint16_t>
class ShortVector 
{
  static_assert(std::is_integral_v<TIndex> && std::is_unsigned_v<TIndex>,
                "TIndex must be an unsigned integral type");
  static_assert(std::is_trivially_copyable_v<TValue> &&
                std::is_trivially_destructible_v<TValue>,
                "TValue must be trivially copyable and destructible");

public:

  using value_type = TValue;
  using size_type  = TIndex;

  ShortVector() noexcept = default;

  explicit ShortVector(size_type count, const value_type& init = value_type{}) 
  {
    reserve(count);
    std::fill_n(data(), count, init);
    setSize_(count);
  }

  ~ShortVector() 
  { 
    deallocate_(); 
  }

  ShortVector(const ShortVector&)            = delete;
  ShortVector& operator=(const ShortVector&) = delete;

  ShortVector(ShortVector&& other) noexcept
  : m_Block(other.m_Block) 
  {
    other.m_Block = nullptr;
  }

  ShortVector& operator=(ShortVector&& other) noexcept 
  {
    if (this != &other) {
      deallocate_();
      m_Block = other.m_Block;
      other.m_Block = nullptr;
    }
    return *this;
  }

  // --- element access (no hidden mutations) --------------------------------
  value_type& operator[](size_type i)             { return data()[i]; }
  const value_type& operator[](size_type i) const { return data()[i]; }

  value_type& at(size_type i) 
  {
    if (i >= size()) throw std::out_of_range("ShortVector::at");
    return (*this)[i];
  }

  const value_type& at(size_type i) const 
  {
    if (i >= size()) throw std::out_of_range("ShortVector::at");
    return (*this)[i];
  }

  value_type*       data()       { return payload_(m_Block); }
  const value_type* data() const { return payload_(m_Block); }

  value_type*       begin()       { return data(); }
  value_type*       end()         { return data() + size(); }
  const value_type* begin() const { return data(); }
  const value_type* end()   const { return data() + size(); }

  // --- capacity -------------------------------------------------------------
  size_type size() const     { return m_Block ? *sizePtr_(m_Block) : size_type(0); }
  size_type capacity() const { return m_Block ? *capPtr_(m_Block)  : size_type(0); }
  bool      empty()   const  { return size() == 0; }

  void reserve(size_type new_cap) {
    if (new_cap <= capacity()) return;
    if (new_cap == 0) return;
    reallocateCopy_(new_cap);
  }

  void shrink_to_fit() {
    if (!m_Block) return;
    const size_type sz = size();
    if (sz == capacity()) return;
    if (sz == 0) { deallocate_(); return; }
    reallocateCopy_(sz);
  }

  void clear() noexcept {
    if (m_Block) setSize_(0);
  }

  // --- modifiers ------------------------------------------------------------
  void resize(size_type new_size, const value_type& init = value_type{}) {
    const size_type old = size();
    if (new_size <= old) { setSize_(new_size); return; }

    ensureGrow_(new_size);
    std::fill(data() + old, data() + new_size, init);
    setSize_(new_size);
  }

  void push_back(const value_type& v) {
    if (size() == capacity()) {
      const size_type cur = capacity();
      const size_type maxN = std::numeric_limits<size_type>::max();
      if (cur == maxN) throw std::length_error("ShortVector capacity limit");
      size_type new_cap = (cur == 0) ? size_type(8) : size_type(std::min<std::uint64_t>(maxN, std::uint64_t(cur) * 2));
      ensureGrow_(new_cap);
    }
    data()[size()] = v;
    setSize_(size() + 1);
  }

  // approximate memory usage: object pointer + header + payload bytes
  std::size_t memoryUsageBytes() const {
    if (!m_Block) return sizeof(m_Block);
    return sizeof(m_Block) + headerBytes_() + std::size_t(capacity()) * sizeof(value_type);
  }

private:

  // Pointer to start of allocated block (header at offset 0)
  std::uint8_t* m_Block = nullptr;

  // ---- header/payload math -------------------------------------------------
  static constexpr std::size_t headerRaw_() noexcept {
    return 2 * sizeof(size_type);
  }
  static constexpr std::size_t headerBytes_() noexcept {
    return alignUp(headerRaw_(), alignOf<value_type>());
  }

  static size_type* sizePtr_(std::uint8_t* blk) {
    return reinterpret_cast<size_type*>(blk + 0);
  }
  static const size_type* sizePtr_(const std::uint8_t* blk) {
    return reinterpret_cast<const size_type*>(blk + 0);
  }

  static size_type* capPtr_(std::uint8_t* blk) {
    return reinterpret_cast<size_type*>(blk + sizeof(size_type));
  }
  static const size_type* capPtr_(const std::uint8_t* blk) {
    return reinterpret_cast<const size_type*>(blk + sizeof(size_type));
  }

  static value_type* payload_(std::uint8_t* blk) {
    return blk ? reinterpret_cast<value_type*>(blk + headerBytes_()) : nullptr;
  }
  static const value_type* payload_(const std::uint8_t* blk) {
    return blk ? reinterpret_cast<const value_type*>(blk + headerBytes_()) : nullptr;
  }

  static std::size_t bytesFor_(size_type cap) {
    return headerBytes_() + std::size_t(cap) * sizeof(value_type);
  }

  // ---- allocation paths (over-aligned if needed) ---------------------------
  static std::uint8_t* alloc_(std::size_t total_bytes) {
    if constexpr (alignOf<value_type>() <= alignof(std::max_align_t)) {
      return static_cast<std::uint8_t*>(::operator new(total_bytes));
    } else {
      return static_cast<std::uint8_t*>(::operator new(total_bytes, std::align_val_t(alignOf<value_type>())));
    }
  }

  static void dealloc_(void* ptr) {
    if (!ptr) return;
    if constexpr (alignOf<value_type>() <= alignof(std::max_align_t)) {
      ::operator delete(ptr);
    } else {
      ::operator delete(ptr, std::align_val_t(alignOf<value_type>()));
    }
  }

  void deallocate_() {
    dealloc_(m_Block);
    m_Block = nullptr;
  }

  void setSize_(size_type n) {
    // caller guarantees 0 <= n <= capacity()
    *sizePtr_(m_Block) = n;
  }

  void ensureGrow_(size_type min_cap) {
    if (min_cap <= capacity()) return;
    reallocateCopy_(min_cap);
  }

  void reallocateCopy_(size_type new_cap) {
    // clamp to max index
    const size_type maxN = std::numeric_limits<size_type>::max();
    if (new_cap > maxN) throw std::length_error("ShortVector capacity limit");

    const std::size_t total = bytesFor_(new_cap);
    std::uint8_t* new_block = alloc_(total);

    // write header
    *reinterpret_cast<size_type*>(new_block + 0) = size();      // size stays same initially
    *reinterpret_cast<size_type*>(new_block + sizeof(size_type)) = new_cap;

    // copy payload if any
    const size_type sz = size();
    if (m_Block && sz) {
      std::memcpy(payload_(new_block), payload_(m_Block), std::size_t(sz) * sizeof(value_type));
    }

    deallocate_();
    m_Block = new_block;
  }
};

} // namespace

// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

TEST_CASE("ShortVector_u32_basics[sv][u32]") {
  edl::ShortVector<std::uint32_t> v;
  CHECK(v.size() == 0);
  CHECK(v.capacity() == 0);
  CHECK(v.empty());
  CHECK(v.data() == nullptr);

  v.push_back(1);
  v.push_back(2);
  v.push_back(3);
  CHECK(v.size() == 3);
  CHECK(v[0] == 1);
  CHECK(v[1] == 2);
  CHECK(v[2] == 3);

  std::uint64_t sum = 0;
  for (auto x : v) sum += x;
  CHECK(sum == 6);

  CHECK_THROWS_AS(v.at(3), std::out_of_range);
  CHECK(v.at(1) == 2);
}

TEST_CASE("ShortVector_u32_reserve_resize_shrink[sv][u32]") {
  edl::ShortVector<std::uint32_t> v;
  v.reserve(10);
  CHECK(v.capacity() >= 10);
  CHECK(v.size() == 0);

  v.resize(5, 7);
  CHECK(v.size() == 5);
  for (auto x : v) CHECK(x == 7);

  const auto cap_before = v.capacity();
  v.shrink_to_fit();
  CHECK(v.capacity() == v.size());
  CHECK(v.capacity() <= cap_before);

  v.clear();
  CHECK(v.size() == 0);
  CHECK(v.capacity() == v.capacity());
}

TEST_CASE("ShortVector_growth_and_max_size_clamp_u8[sv][u32][idx8]") {
  using SV = edl::ShortVector<std::uint32_t, std::uint8_t>;
  SV v;
  for (int i = 0; i < 200; ++i) v.push_back(static_cast<std::uint32_t>(i));
  CHECK(v.size() == 200);
  CHECK(v.capacity() >= 200);
  CHECK(v[199] == 199);

  v.reserve(std::numeric_limits<std::uint8_t>::max());
  CHECK(v.capacity() == std::numeric_limits<std::uint8_t>::max());

  SV w;
  w.reserve(std::numeric_limits<std::uint8_t>::max());
  w.resize(std::numeric_limits<std::uint8_t>::max());
  CHECK_THROWS_AS(w.push_back(42), std::length_error);
}

TEST_CASE("ShortVector_u64_alignment[sv][u64]") {
  edl::ShortVector<std::uint64_t> v;
  for (int i = 0; i < 9; ++i) v.push_back(0xDEADBEEFCAFEBABEull + i);
  CHECK(v.size() == 9);
  CHECK(reinterpret_cast<std::uintptr_t>(v.data()) % alignof(std::uint64_t) == 0);
  CHECK(v[0] == 0xDEADBEEFCAFEBABEull);
  CHECK(v[8] == 0xDEADBEEFCAFEBABEull + 8);
}

TEST_CASE("ShortVector_move_only_semantics[sv][u32][move]") {
  static_assert(!std::is_copy_constructible_v<edl::ShortVector<std::uint32_t>>);
  static_assert(!std::is_copy_assignable_v<edl::ShortVector<std::uint32_t>>);

  edl::ShortVector<std::uint32_t> a;
  for (int i = 0; i < 16; ++i) a.push_back(i);

  edl::ShortVector<std::uint32_t> b(std::move(a));
  CHECK(a.size() == 0);
  CHECK(b.size() == 16);
  CHECK(b[15] == 15);

  edl::ShortVector<std::uint32_t> c;
  c = std::move(b);
  CHECK(b.size() == 0);
  CHECK(c.size() == 16);
  CHECK(c[0] == 0);
}

TEST_CASE("ShortVector_u32_u16_memoryUsageBytes[sv][u32][idx16]") {
  edl::ShortVector<std::uint32_t> v;
  auto base = v.memoryUsageBytes();
  v.reserve(10);
  auto after = v.memoryUsageBytes();
  CHECK(after >= base + sizeof(std::uint32_t) * 10);
}

#include <random>

TEST_CASE("ShortVector_u32_u16_item_loop[sv][u32][idx16]") 
{
  using namespace EDL_NAMESPACE;
  using namespace std;
  //
  typedef ShortVector<uint32_t,uint16_t> shortvec_t;
  const size_t N = 10000;
  //
  // create a std::vector with N random values of 32 bit integers
  //
  vector<uint32_t> data(N);
  random_device rd;
  mt19937_64 eng(rd());
  uniform_int_distribution<uint32_t> distr(0, numeric_limits<uint32_t>::max());
  for (size_t i = 0; i < N; ++i) {
    data[i] = distr(eng);
  }    
  //
  // copy the data to a ShortVector
  //
  shortvec_t sv;
  for (size_t i = 0; i < N; ++i) {
    sv.push_back(data[i]);
  }
  //
  // check if the data is correct (index access)
  //
  int i = 0;
  for (auto v : sv) {
    CHECK(v == data[i]);
    ++i;
  }
}
