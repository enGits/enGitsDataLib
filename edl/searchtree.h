// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef SEARCHTREE_H
#define SEARCHTREE_H

#include <cstdint>
#include <type_traits>
#include <limits>
#include <vector>

#include "edl/amr.h"
#include "edl/edl.h"
#include "edl/mathvector.h"


namespace edl
{

template <typename VEC> struct SearchTreeVectorCheck 
{
  typedef VEC item_type;
  static bool isInsideCartesianBox(const VEC& x, const VEC& x1, const VEC& x2)
  {
    for (int dim = 0; dim < x.size(); ++dim) {
      if (x[dim] <= x1[dim] || x[dim] >= x2[dim]) {
        return false;
      }
    }
    return true;
  }
  static bool match(const VEC& x1, const VEC& x2)
  {
    return false;
  }
  static void boundingBox(const VEC& item, VEC& x1, VEC& x2)
  {
    x1 = item;
    x2 = item;
  }
};


template <typename TIndex, typename Tijk, typename TVector, typename TCheck>
class TSearchTree
{
  public: // data types

    typedef          TIndex             index_type;
    typedef          TVector            vector_type;
    typedef          Tijk               ijk_type;
    typedef          AMRIndex<ijk_type> amr_index_type;
    typedef          TCheck             check_type;
    typedef typename TCheck::item_type  item_type;


protected: // attributes

  AMRMesh<TIndex, Tijk, TVector>*                      m_AMR = nullptr;
  vector_type                                          m_X1;
  vector_type                                          m_X2;
  real                                                 m_H;
  int                                                  m_MaxLevel = 10;
  int                                                  m_MaxBucketSize;
  std::unordered_map<amr_index_type, std::vector<int>> m_Buckets;
  std::vector<item_type>                               m_Items;
  real                                                 m_MaxSearchDist = small_real;

  static const int                                     m_MaxAllowedLevel = 10;


public: // methods

  TSearchTree(int max_bucket_size=10) : m_MaxBucketSize(max_bucket_size) {};

  ~TSearchTree()
  {
    if (m_AMR) {
      delete m_AMR;
    }
  }

  void setMaximalBucketSize(int max_bucket_size)
  {
    m_MaxBucketSize = max_bucket_size;
  }

  void setMaxSearchDist(real max_search_dist)
  {
    m_MaxSearchDist = max_search_dist;
  }

  int computeMaxLevel(real H)
  {
    using namespace std;
    //
    auto n = H/m_MaxSearchDist;
    int  max_level = 0;
    while (max_level < m_MaxAllowedLevel) {
      auto l2 = ipow(2, max_level+1);
      if (l2 > n) {
        break;
      }
      ++max_level;
    }
    return max_level;
  }

  template <typename C>
  void setItems(const C items)
  {
    //
    // check that C::value_type is item_type
    //
    static_assert(std::is_same<typename C::value_type, item_type>::value, "C::value_type must be item_type");
    //
    // compute the bounding box of all items
    //
    m_X1 = vector_type( std::numeric_limits<typename vector_type::value_type>::max());
    m_X2 = vector_type(-std::numeric_limits<typename vector_type::value_type>::max());
    for (auto item : items) {
      vector_type x1, x2;
      check_type::boundingBox(item, x1, x2);
      for (int i = 0; i < 3; ++i) {
        m_X1[i] = std::min(m_X1[i], x1[i]);
        m_X2[i] = std::max(m_X2[i], x2[i]);
      }
    }
    auto xc = 0.5 * (m_X1 + m_X2);
    auto d  = 0.5 * (m_X2 - m_X1);
    m_H  = std::max(d[0], std::max(d[1], d[2]));
    m_X1 = xc - vector_type(m_H);
    m_X2 = xc + vector_type(m_H);
    //
    vector_type delta = (m_X2 - m_X1);
    delta.normalise();
    static real eps = 1e-6;
    m_X1 -= eps*m_H + m_MaxSearchDist;
    m_X2 += eps*m_H + m_MaxSearchDist;
    //
    // set up the AMR mesh
    //
    m_Items = items;
    m_AMR   = new AMRMesh<TIndex, Tijk, TVector>(1, 1, 1, m_X1, m_X2);
    //
    amr_index_type root(0,0,0,0);
    m_Buckets[root] = std::vector<int>(m_Items.size());
    for (int i = 0; i < m_Items.size(); ++i) {
      m_Buckets[root][i] = i;
    }
    m_H  = std::max(m_X2[0] - m_X1[0], std::max(m_X2[1] - m_X1[1], m_X2[2] - m_X1[2]));
    //
    // recursively refine the AMR mesh
    //
    int max_level = 0;
    m_MaxLevel = std::min(computeMaxLevel(m_H), m_MaxLevel);
    if (m_Items.size() > m_MaxBucketSize) {
      m_AMR->refineCell(root);
      for (int level = 1; level <= m_MaxLevel; ++level) {
        auto leaf_cells = m_AMR->getLeafCellIndices();
        bool refined    = false;
        for (auto ijk : leaf_cells) {
          if (ijk.level() == level) {
            auto h  = m_H/ipow(2, level);
            auto bb = m_AMR->getCellBoundingBox(ijk);
            auto x1 = bb.first;
            auto x2 = bb.second;
            x1 -= eps*h + m_MaxSearchDist; 
            x2 += eps*h + m_MaxSearchDist;
            for (auto i : m_Buckets[ijk.parent()]) {
              if (check_type::isInsideCartesianBox(m_Items[i], x1, x2)) {
                m_Buckets[ijk].push_back(i);
              }
            }
            if (m_Buckets[ijk].size() > m_MaxBucketSize) {
              if (level < m_MaxLevel) {
                m_AMR->refineCell(ijk);
                refined = true;
              }
            }
          }
        }
        if (!refined) {
          break;
        }
        max_level = level;
      }
    }
    //
    // clean up (i.e. delete all buckets which are not leaf cells)
    //
    int N = 1;
    for (int level = 0; level <= max_level; ++level) {
      for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
          for (int k = 0; k < N; ++k) {
            amr_index_type ijk(i,j,k,level);
            if (!m_AMR->isLeafCell(ijk)) {
              if (m_Buckets.find(ijk) != m_Buckets.end()) {
                m_Buckets.erase(ijk);
              }
            }
          }
        }
      }
      N *= 2;
    }
  }

  amr_index_type getCellContainingPointOnLevel(const vector_type& x, int level) const
  {
    using namespace std;
    auto h = m_H/ipow(2, level);
    int  i = (x[0] - m_X1[0])/h;
    int  j = (x[1] - m_X1[1])/h;
    int  k = (x[2] - m_X1[2])/h;
    return amr_index_type(i, j, k, level);
  }

  int nearestItemIndex(const vector_type& x) const
  {
    int  idx = -1;
    auto ijk = amr_index_type(0,0,0,0);
    if (check_type::isInsideCartesianBox(x, m_X1, m_X2)) {
      while (!m_AMR->isLeafCell(ijk)) {
        ijk = getCellContainingPointOnLevel(x, ijk.level() + 1);
      }
    }
    real min_dist = max_real;
    if (m_Buckets.find(ijk) != m_Buckets.end()) {
      for (auto i : m_Buckets.at(ijk)) {
        vector_type x1, x2;
        check_type::boundingBox(m_Items[i], x1, x2);
        auto xc = 0.5*(x1 + x2);
        auto d  = (x - xc).abs();
        if (d < min_dist) {
          min_dist = d;
          idx = i;
        }
      }
    }
    return idx;
  }

  int findMatchingItem(const vector_type& x) const
  {
    return -1;
  }

};

template <class T>
using PointSearchTree = TSearchTree<uint32_t, uint16_t, MathVector<StaticVector<T,3>>, SearchTreeVectorCheck<MathVector<StaticVector<T,3>>>>;

} // namespace edl

// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

#include <random>
#include "edl/quicktimer.h"
#include "edl/geometrytools.h"

TEST_CASE("SearchTree__compute_max_level")
{
  using namespace edl;
  //
  PointSearchTree<real> tree(10);
  tree.setMaxSearchDist(1.0);
  CHECK(tree.computeMaxLevel(1.00)  == 0);
  CHECK(tree.computeMaxLevel(1.99)  == 0);
  CHECK(tree.computeMaxLevel(2.00)  == 1);
  CHECK(tree.computeMaxLevel(2.01)  == 1);
  CHECK(tree.computeMaxLevel(3.00)  == 1);
  CHECK(tree.computeMaxLevel(4.00)  == 2);
  CHECK(tree.computeMaxLevel(5.00)  == 2);
  CHECK(tree.computeMaxLevel(6.00)  == 2);
  CHECK(tree.computeMaxLevel(7.00)  == 2);
  CHECK(tree.computeMaxLevel(8.00)  == 3);
  CHECK(tree.computeMaxLevel(9.00)  == 3);
  CHECK(tree.computeMaxLevel(10.00) == 3);
  CHECK(tree.computeMaxLevel(11.00) == 3);
  CHECK(tree.computeMaxLevel(12.00) == 3);
  CHECK(tree.computeMaxLevel(13.00) == 3);
  CHECK(tree.computeMaxLevel(14.00) == 3);
  CHECK(tree.computeMaxLevel(15.00) == 3);
  CHECK(tree.computeMaxLevel(16.00) == 4);
  CHECK(tree.computeMaxLevel(17.00) == 4);
}


TEST_CASE("SearchTree__random_Items_search")
{
  using namespace std;
  typedef edl::real real;
  typedef edl::MathVector<edl::StaticVector<real,3>> vec_t;
  //
  const int  num_items = 2000;
  const int  num_tests = 100;
  const real range     = 3;
  QuickTimer timer;
  //
  // create a set of random items on a sphere with radius range/3
  //
  mt19937 gen(42);
  uniform_real_distribution<real> dist1(-M_PI, M_PI);
  vector<vec_t> items(num_items);
  vec_t y_axis(0,1,0);
  vec_t z_axis(0,0,1);
  for (int i = 0; i < num_items; ++i) {
    items[i] = vec_t(range/3,0,0);
    items[i] = rotate(items[i], y_axis, dist1(gen));
    items[i] = rotate(items[i], z_axis, dist1(gen));
  }
  //
  // create the search tree
  //
  timer.start();
  edl::PointSearchTree<real> search(10);
  real max_search_dist = range/10;
  search.setMaxSearchDist(max_search_dist);
  search.setItems(items);
  timer.stop();
  cout << "SearchTree construction time : " << timer.milliseconds() << " ms" << endl;
  //
  // create a set of test items
  //
  uniform_real_distribution<real> dist2(-range/2, range/2);
  vector<vec_t> test_items(num_tests);
  for (int i = 0; i < num_tests; ++i) {
    test_items[i] = vec_t(dist2(gen), dist2(gen), dist2(gen));
  }
  test_items.push_back(vec_t(-range/2,-range/2,-range/2));
  test_items.push_back(vec_t( range/2, range/2, range/2));
  //
  // brute force loop to find the reference solution
  //
  timer.restart();
  vector<int> nearest_item(num_tests);
  int N_found = 0;  
  for (int i = 0; i < test_items.size(); ++i) {
    vec_t p = test_items[i];
    real  d_min = max_search_dist;
    nearest_item[i] = -1;
    for (int j = 0; j < num_items; ++j) {
      vec_t q = items[j];
      real  d = (p - q).abs();
      if (d < d_min) {
        d_min = d;
        nearest_item[i] = j;        
      }          
    }
    if (nearest_item[i] >= 0) {
      ++N_found;
    }
  }
  timer.stop();
  cout << "Brute force search time      : " << timer.milliseconds() << " ms" << endl;
  cout << "Number of found items        : " << N_found << " of " << num_tests << endl;
  //
  // perform the tests
  //
  timer.restart();
  for (int i = 0; i < test_items.size(); ++i) {
    vec_t p = test_items[i];
    int   j = search.nearestItemIndex(p);
    if (nearest_item[i] >= 0) {
      auto d = (items[nearest_item[i]] - p).abs();
      CHECK(j == nearest_item[i]);
    }
  }
  timer.stop();
  cout << "SearchTree search time       : " << timer.milliseconds() << " ms" << endl;
}


#endif // SEARCHTREE_H