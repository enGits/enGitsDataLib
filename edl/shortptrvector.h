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
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include "shortvector.h"

namespace EDL_NAMESPACE
{

// replace the previous detection helpers with these C++17 versions:

template <class, class = void>
struct has_data : std::false_type {};
template <class C>
struct has_data<C, std::void_t<decltype(std::declval<const C&>().data())>> : std::true_type {};
template <class C>
inline constexpr bool has_data_v = has_data<C>::value;

template <class, class = void>
struct has_value_type : std::false_type {};
template <class C>
struct has_value_type<C, std::void_t<typename C::value_type>> : std::true_type {};
template <class C>
inline constexpr bool has_value_type_v = has_value_type<C>::value;

// Encodes pointers as element-offsets into a single, shared backing container.
// Per-instance storage == ShortVector<offsets>. The backing container pointer is static.
template<
  class T,
  class TIndex = std::uint16_t,
  class TOff   = std::uint32_t,
  class TContainer = std::vector<T>
>
class ShortPtrVector {
  static_assert(std::is_unsigned_v<TIndex> && std::is_integral_v<TIndex>, "TIndex must be an unsigned integral type");
  static_assert(std::is_unsigned_v<TOff>   && std::is_integral_v<TOff>,   "TOff must be an unsigned integral type");
  static_assert(!has_value_type_v<TContainer> || std::is_same_v<typename TContainer::value_type, T>, "TContainer::value_type must be T");
public:
  using value_type       = T*;
  using size_type        = TIndex;
  using off_type         = TOff;
  using difference_type  = std::ptrdiff_t;
  using container_type   = TContainer;

private:
  static constexpr off_type kNullOff = std::numeric_limits<off_type>::max();

  // Shared (non-owning) pointer to the backing container. Must outlive all ShortPtrVector instances.
  inline static const container_type* m_Container = nullptr;

  // The only per-instance state: compact offsets
  ShortVector<off_type, TIndex> m_Offsets;

  static const T* basePtr_() {
    if (!m_Container) return nullptr;
    if constexpr (has_data_v<container_type>) {
      return m_Container->size() ? m_Container->data() : nullptr;
    } else {
      return m_Container->size() ? std::addressof((*m_Container)[0]) : nullptr;
    }
  }

  static off_type encode_(value_type p) {
    if (p == nullptr) return kNullOff;
    const T* base = basePtr_();
    if (!base) throw std::logic_error("ShortPtrVector: container not set or empty");
    const auto diff = p - base; // elements
    if (diff < 0) throw std::out_of_range("ShortPtrVector: pointer before container base");
    using UDiff = std::make_unsigned_t<difference_type>;
    const auto udiff = static_cast<UDiff>(diff);
    if (udiff > std::numeric_limits<off_type>::max())
      throw std::out_of_range("ShortPtrVector: offset does not fit in off_type");
    return static_cast<off_type>(udiff);
  }

  static value_type decode_(off_type off) {
    if (off == kNullOff) return nullptr;
    if (!m_Container) throw std::logic_error("ShortPtrVector: container not set");
    // No need for base: contiguous + operator[] gives us a reference; take its address.
    return const_cast<T*>(std::addressof((*m_Container)[static_cast<std::size_t>(off)]));
  }

public:
  // --- lifecycle (move-only) ------------------------------------------------
  ShortPtrVector() = default;
  explicit ShortPtrVector(size_type count, value_type init = nullptr) { resize(count, init); }

  ShortPtrVector(const ShortPtrVector&) = delete;
  ShortPtrVector& operator=(const ShortPtrVector&) = delete;

  ShortPtrVector(ShortPtrVector&& other) noexcept
  : m_Offsets(std::move(other.m_Offsets)) {}

  ShortPtrVector& operator=(ShortPtrVector&& other) noexcept {
    if (this != &other) m_Offsets = std::move(other.m_Offsets);
    return *this;
  }

  // --- container binding (static) -------------------------------------------
  static void setContainer(const container_type& c) noexcept { m_Container = &c; }
  static const container_type* container() noexcept { return m_Container; }

  // --- capacity -------------------------------------------------------------
  size_type size()     const noexcept { return m_Offsets.size(); }
  size_type capacity() const noexcept { return m_Offsets.capacity(); }
  bool      empty()    const noexcept { return m_Offsets.empty(); }

  void reserve(size_type n)      { m_Offsets.reserve(n); }
  void shrink_to_fit()           { m_Offsets.shrink_to_fit(); }
  void clear() noexcept          { m_Offsets.clear(); }

  // --- element access (return BY VALUE) -------------------------------------
  value_type operator[](size_type i) const { return decode_(m_Offsets[i]); }
  value_type at(size_type i) const         { return decode_(m_Offsets.at(i)); }

  // --- modifiers ------------------------------------------------------------
  void push_back(value_type p) { m_Offsets.push_back(encode_(p)); }

  // Convenience: push by index into the shared container
  void push_back_index(std::size_t idx) {
    // Optional: bounds assert on container size
    m_Offsets.push_back(static_cast<off_type>(idx));
  }

  void resize(size_type n, value_type init = nullptr) {
    const auto old = m_Offsets.size();
    if (n <= old) { m_Offsets.resize(n); return; }
    m_Offsets.resize(n);
    const off_type enc = encode_(init);
    std::fill(m_Offsets.begin() + old, m_Offsets.end(), enc);
  }

  // Debugging: inspect raw offsets
  const ShortVector<off_type, TIndex>& offsets() const noexcept { return m_Offsets; }

  // Approximate memory usage (excludes static container pointer)
  std::size_t memoryUsageBytes() const {
    return sizeof(*this) - sizeof(m_Offsets) + m_Offsets.memoryUsageBytes();
  }
};

} // namespace EDL_NAMESPACE

// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------


// A minimal contiguous container with operator[] and size(), but NO .data()
template<class T>
struct FauxVec {
  using value_type = T;
  std::vector<T> impl;
  std::size_t size() const { return impl.size(); }
  const T& operator[](std::size_t i) const { return impl[i]; }
  T&       operator[](std::size_t i)       { return impl[i]; }
};

TEST_CASE("ShortPtrVector_std_vector_realloc_ok[idx16][off32]") 
{
  using namespace EDL_NAMESPACE;

  using PV = ShortPtrVector<std::uint32_t, std::uint16_t, std::uint32_t, std::vector<std::uint32_t>>;
  std::vector<std::uint32_t> backing;
  backing.reserve(4);
  for (int i=0;i<4;++i) backing.push_back(i);
  PV::setContainer(backing);

  PV pv;
  pv.push_back(&backing[0]);
  pv.push_back(&backing[3]);

  // force reallocation
  for (int i=0;i<2000;++i) backing.push_back(i);

  CHECK(pv[0] == &backing[0]);
  CHECK(pv[1] == &backing[3]);
}

TEST_CASE("ShortPtrVector_custom_container_no_data[idx16][off32]") 
{
  using namespace EDL_NAMESPACE;

  using PV = ShortPtrVector<int, std::uint16_t, std::uint32_t, FauxVec<int>>;
  FauxVec<int> cont;
  cont.impl.resize(10);
  PV::setContainer(cont);

  PV pv;
  pv.push_back(&cont.impl[2]);
  pv.push_back(&cont.impl[9]);
  CHECK(pv[0] == &cont.impl[2]);
  CHECK(pv[1] == &cont.impl[9]);
}

TEST_CASE("ShortPtrVector_nullptr_and_bounds[idx8][off8]") 
{
  using namespace EDL_NAMESPACE;

  using PV = ShortPtrVector<int, std::uint8_t, std::uint8_t, std::vector<int>>;
  std::vector<int> backing(300);
  PV::setContainer(backing);

  PV pv;
  pv.push_back(nullptr);
  CHECK(pv[0] == nullptr);
  pv.push_back(&backing[255]);
  CHECK_THROWS_AS(pv.push_back(&backing[256]), std::out_of_range); // off8 overflow
}

#include <random>

TEST_CASE("ShortPtrVector_write_access[idx16][off32]")
{
  using namespace EDL_NAMESPACE;
  using namespace std;

  int N = 10000;
  vector<uint64_t> data(N);
  random_device rd;
  mt19937_64 eng(rd());
  uniform_int_distribution<uint64_t> distr(0, numeric_limits<uint64_t>::max());
  for (size_t i = 0; i < N; ++i) {
    data[i] = distr(eng);
  }    

  using PV = ShortPtrVector<uint64_t>;
  PV::setContainer(data);
  PV pv;
  for (int i = 0; i < N; ++i) {
    pv.push_back(&data[i]);
  }

  for (int i = 0; i < N; ++i) {
    CHECK(&data[i] == pv[i]);
  }
}

