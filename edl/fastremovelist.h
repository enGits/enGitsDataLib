// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef FASTREMOVELIST_H
#define FASTREMOVELIST_H

#include "edl.h"
#include <vector>

#include <iostream>

namespace EDL_NAMESPACE
{
  
template <typename C>
class FastRemoveList
{
  typedef typename C::value_type value_type;

  const C& m_Container;

  struct node_t
  {
    int prev;
    int next;
  };

  int m_First;
  int m_Last;

  std::vector<node_t> m_Nodes;

public:

  FastRemoveList(const C& container) : m_Container(container)
  {
    m_Nodes.resize(container.size());
    for (int i = 0; i < container.size(); ++i) {
      m_Nodes[i].prev = i - 1;
      m_Nodes[i].next = i + 1;
    }
    m_First = 0;
    m_Last = container.size() - 1;
  }

  void remove(int i)
  {
    int prev = m_Nodes[i].prev;
    int next = m_Nodes[i].next;
    if (prev == -1) {
      m_First = next;
    } else {
      m_Nodes[prev].next = next;
    }
    if (next >= m_Nodes.size()) {
      m_Last = prev;
    } else {
      m_Nodes[next].prev = prev;
    }
    m_Nodes[i].prev = -1;
    m_Nodes[i].next = -1;
  }

  bool hasBeenRemoved(int i) const { return m_Nodes[i].prev == -1 && m_Nodes[i].next == -1; }
  int  firstIndex()          const { return m_First; }
  int  nextIndex(int i)      const { return m_Nodes[i].next; }
  int  lastIndex()           const { return m_Last; }

  const value_type& operator[](int i) const { return m_Container[i]; }
  const value_type& first()           const { return m_Container[m_First]; }
  const value_type& last()            const { return m_Container[m_Last]; }

  value_type popFirst()
  {
    value_type result = first();
    remove(m_First);
    return result;
  }

  bool isEmpty() const { return m_First == m_Container.size(); }

};

} // namespace EDL_NAMESPACE

TEST_CASE("FastRemoveList 01")
{
  #include <vector>
  int N = 10;
  std::vector<int> v(10);
  edl::FastRemoveList<std::vector<int>> frl(v);
  CHECK(frl.firstIndex() == 0);
  CHECK(frl.lastIndex() == N - 1);
  for (int i = 0; i < N; ++i) {
    CHECK(frl.nextIndex(i) == i + 1);
  }
  for (int i = 0; i < N-1; ++i) {
    CHECK(!frl.hasBeenRemoved(i));
    frl.remove(i);
    CHECK(!frl.isEmpty());
    CHECK(frl.hasBeenRemoved(i));
    CHECK(frl.firstIndex() == i + 1);
    CHECK(frl.lastIndex() == N - 1);    
  }
  frl.remove(N-1);
  CHECK(frl.firstIndex() == N);
  CHECK(frl.lastIndex() == -1);
  CHECK(frl.isEmpty());
}

#endif // FASTREMOVELIST_H