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
#include <fstream>
#include <type_traits>
#include <limits>
#include <vector>
#include <filesystem>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <string>

#include "edl/amr.h"
#include "edl/edl.h"
#include "edl/mathvector.h"
#include "edl/geometrytools.h"


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

  void setMaxBucketSize(int max_bucket_size)
  {
    m_MaxBucketSize = max_bucket_size;
  }

  void setMaxSearchDist(real max_search_dist)
  {
    m_MaxSearchDist = max_search_dist;
  }

  int computeMaxLevel(real H)
  {
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
  // void setItems(const C items)
  void setItems(const C& items)
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
    for (item_type item : items) {
      vector_type x1, x2;
      check_type::boundingBox(item, x1, x2);
      for (int i = 0; i < 3; ++i) {
        m_X1[i] = std::min(m_X1[i], x1[i]);
        m_X2[i] = std::max(m_X2[i], x2[i]);
      }
    }
    vector_type xc = 0.5 * (m_X1 + m_X2);
    vector_type d  = 0.5 * (m_X2 - m_X1);
    m_H  = std::max(d[0], std::max(d[1], d[2]));
    m_X1 = xc - vector_type(m_H);
    m_X2 = xc + vector_type(m_H);
    //
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
        for (amr_index_type ijk : leaf_cells) {
          if (ijk.level() == level) {
            real h  = m_H/ipow(2, level);
            auto bb = m_AMR->getCellBoundingBox(ijk);
            vector_type x1 = bb.first;
            vector_type x2 = bb.second;
            x1 -= eps*h + m_MaxSearchDist; 
            x2 += eps*h + m_MaxSearchDist;
            for (int i : m_Buckets[ijk.parent()]) {
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
    auto        h = m_H/ipow(2, level);
    // Unused locals ? "x0, x1, x2" (fatal under -Werror) 
    vector_type x0 = (x[0] - m_X1[0])/h;
    vector_type x1 = (x[1] - m_X1[1])/h;
    vector_type x2 = (x[2] - m_X1[2])/h;

    int         i  = (x[0] - m_X1[0])/h;
    int         j  = (x[1] - m_X1[1])/h;
    int         k  = (x[2] - m_X1[2])/h;
    return amr_index_type(i, j, k, level);
  }

  std::vector<int> getBucketContents(const vector_type& x) const
  {
    std::vector<int> items;
    amr_index_type ijk(0,0,0,0);
    if (check_type::isInsideCartesianBox(x, m_X1, m_X2)) {
      while (!m_AMR->isLeafCell(ijk)) {
        ijk = getCellContainingPointOnLevel(x, ijk.level() + 1);
      }
    }
    if (m_Buckets.find(ijk) != m_Buckets.end()) {
      items = m_Buckets.at(ijk);
    }
    return items;
  }

  int nearestItemIndex(const vector_type& x) const
  {
    int  idx = -1;
    real min_dist = max_real;
    for (int i : getBucketContents(x)) {
      vector_type x1, x2;
      check_type::boundingBox(m_Items[i], x1, x2);
      vector_type xc = 0.5*(x1 + x2);
      real        d  = (x - xc).abs();
      if (d < min_dist) {
        min_dist = d;
        idx = i;
      }
    }
    return idx;
  }

  int findMatchingItem(const vector_type& x, std::vector<int> items = {}) const
  {
    if (items.empty()) {
      items = getBucketContents(x);
    }
    for (int item_index : items) {
      if (check_type::match(x, m_Items[item_index])) {
        return item_index;
      }
    }
    return -1;
  }

  const item_type& getItem(int index) const
  {
    return m_Items[index];
  }

  void writeFoamDebugMesh(std::string path)
  {
    //
    // make sure the path exists
    //    
    namespace fs = std::filesystem;
    fs::path root_dir  = path;
    fs::path const_dir = root_dir / "constant";
    fs::path poly_dir  = const_dir / "polyMesh";
    if (!fs::exists(root_dir)) {
      fs::create_directories(root_dir);
    }
    if (!fs::exists(const_dir)) {
      fs::create_directories(const_dir);      
    }
    if (!fs::exists(poly_dir)) {
      fs::create_directories(poly_dir);      
    }
    //
    // write the mesh
    //
    m_AMR->writeFoamMesh(poly_dir.string());
    //
    // create .foam file
    std::ofstream ofs(path + "/search_tree.foam");
    ofs << "\n" << std::endl;
  }

  void writeVtkDebugMesh(std::string file_name)
  {
    using std::vector;
    //
    auto leaf_cells = m_AMR->getLeafCellIndices();
    auto points     = m_AMR->extractPoints();
    //
    vector<vector<int>> cells;
    cells.reserve(leaf_cells.size());
    //
    for (auto ijk : leaf_cells) {
      vector<amr_index_type> nodes;
      nodes.reserve(8);
      nodes.push_back(ijk);
      nodes.push_back(ijk.incrementedI());
      nodes.push_back(ijk.incrementedI().incrementedJ());
      nodes.push_back(ijk.incrementedJ());
      nodes.push_back(ijk.incrementedK());
      nodes.push_back(ijk.incrementedI().incrementedK());
      nodes.push_back(ijk.incrementedI().incrementedJ().incrementedK());
      nodes.push_back(ijk.incrementedJ().incrementedK());
      vector<int> cell;
      cell.reserve(8);
      for (auto node : nodes) {
        cell.push_back(m_AMR->nodeLinearIndex(node));
      }
      cells.push_back(cell);
    }
    //
    std::ofstream ofs(file_name);
    if (!ofs) {
      std::cerr << "Error opening file: " << file_name << std::endl;
      return;
    }
  
    // Header
    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "Hex mesh\n";
    ofs << "ASCII\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";
  
    // Points
    ofs << "POINTS " << points.size() << " float\n";
    for (const auto& p : points)
      ofs << p[0] << " " << p[1] << " " << p[2] << "\n";
  
    // Cells
    const int num_cells = cells.size();
    const int entries_per_cell = 9; // 1 count + 8 point indices
    ofs << "CELLS " << num_cells << " " << num_cells * entries_per_cell << "\n";
    for (const auto& cell : cells) {
      if (cell.size() != 8) {
        std::cerr << "Error: Hex cell must have 8 points.\n";
        return;
      }
      ofs << "8 ";
      for (int idx : cell)
        ofs << idx << " ";
      ofs << "\n";
    }
  
    // Cell types (VTK_HEXAHEDRON = 12)
    ofs << "CELL_TYPES " << num_cells << "\n";
    for (int i = 0; i < num_cells; ++i)
      ofs << "12\n";
  
    ofs.close();
    std::cout << "Wrote " << file_name << " successfully.\n";
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
  using std::vector;
  typedef edl::real real;
  typedef edl::MathVector<edl::StaticVector<real,3>> vec_t;
  //
  const int  num_items = 2000;
  const int  num_tests = 100;
  const real range     = 3;
  edl::QuickTimer timer;
  //
  // create a set of random items on a sphere with radius range/3
  //
  std::mt19937 gen(42);
  std::uniform_real_distribution<real> dist1(-M_PI, M_PI);
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
  std::cout << "SearchTree construction time : " << timer.milliseconds() << " ms" << std::endl;
  //
  // create a set of test items
  //
  std::uniform_real_distribution<real> dist2(-range/2, range/2);
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
  std::cout << "Brute force search time      : " << timer.milliseconds() << " ms" << std::endl;
  std::cout << "Number of found items        : " << N_found << " of " << num_tests << std::endl;
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
  std::cout << "SearchTree search time       : " << timer.milliseconds() << " ms" << std::endl;
}


TEST_CASE("SearchTree__points_on_a_line_search")
{
  using std::vector;
  typedef float real;
  typedef edl::MathVector<edl::StaticVector<real,3>> vec_t;
  //
  const int  num_items = 64;
  vec_t x1(-0.5,-0.5,0.5);
  vec_t x2( 0.5,-0.5,0.5);
  real search_dist = 2*0.03125;
  //
  vector<vec_t> items(num_items);
  vec_t x_to_snap(0.015625, -0.499222696, 0.499226183);
  real min_dist;
  int  item_index = -1;
  //
  for (int i = 0; i < num_items; ++i) {
    real x = x1[0] + (x2[0] - x1[0])*i/(num_items-1);
    real y = x1[1] + (x2[1] - x1[1])*i/(num_items-1);
    real z = x1[2] + (x2[2] - x1[2])*i/(num_items-1);
    items[i] = vec_t(x,y,z);
    //
    real d = (x_to_snap - items[i]).abs();
    if (d < min_dist || item_index < 0) {
      min_dist   = d;
      item_index = i;
    }
  }
  //
  // create the search tree
  //
  edl::PointSearchTree<real> search(20);
  search.setMaxSearchDist(search_dist);
  search.setItems(items);
  search.writeVtkDebugMesh("search_tree.vtk");
  //
  int   j = search.nearestItemIndex(x_to_snap);
  CHECK(j == item_index);
}


#endif // SEARCHTREE_H
