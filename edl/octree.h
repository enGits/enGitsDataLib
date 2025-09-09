// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef OCTREE_H
#define OCTREE_H

#include <algorithm>
#include <cmath>
#include <stdint.h>
#include <cstdint>
#include <fstream>
#include <iostream>

#include <ostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>

#include "edl/edl.h"
#include "edl/edlerror.h"
#include "edl/geometrytools.h"
#include "edl/mathvector.h"
#include "edl/quicktimer.h"

namespace edl
{

/*

      6-------------7
     /|            /|
    / |           / |
   /  |          /  |
  /   |         /   |
 /    |        /    |
4-------------5     |
|     |       |     |
|     2-------|-----3
|    /        |    /
|   /         |   /
|  /          |  /
| /           | /
|/            |/
0-------------1

*/

template <typename VEC> struct OctreeVectorCheck 
{
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

template <typename VEC, typename ITEM, typename CHECK>
class Octree
{

public: // data types

  typedef VEC vector_t;
  typedef typename vector_t::value_type scalar_t;
  typedef ITEM item_t;
  typedef CHECK check_t;


protected: // data types

  /**
    * @brief A geometric vertex in the octree.
    */
  struct vertex_t
  {
    std::uint16_t m_Ix = 0;
    std::uint16_t m_Iy = 0;
    std::uint16_t m_Iz = 0;
    std::int32_t  m_Index = -1;
    //
    vertex_t() {};
    //
    vertex_t(vector_t x, const Octree<VEC,ITEM,CHECK>& octree)
    {
      vector_t x1 = octree.getBoundingBoxX1();
      vector_t x2 = octree.getBoundingBoxX2();
      vector_t dx = x2 - x1;
      //
      m_Ix = std::uint16_t(65535.0 * (x[0] - x1[0]) / dx[0]);
      m_Iy = std::uint16_t(65535.0 * (x[1] - x1[1]) / dx[1]);
      m_Iz = std::uint16_t(65535.0 * (x[2] - x1[2]) / dx[2]);
    };
    //
    vertex_t(const vertex_t& v1, const vertex_t& v2)
    {
      m_Ix = (v1.m_Ix + v2.m_Ix) / 2;
      m_Iy = (v1.m_Iy + v2.m_Iy) / 2;
      m_Iz = (v1.m_Iz + v2.m_Iz) / 2;
    };
    //
    vertex_t(uint16_t ix, uint16_t iy, uint16_t iz) : m_Ix(ix), m_Iy(iy), m_Iz(iz) {};
    //
    vector_t getCoordinates(const Octree<VEC,ITEM,CHECK>& octree) const
    {
      scalar_t x1[3], x2[3], dx[3];
      //vector_t x1 = octree.getBoundingBoxX1();
      //vector_t x2 = octree.getBoundingBoxX2();
      //vector_t dx = x2 - x1;
      x1[0] = octree.getBoundingBoxX1()[0];
      x1[1] = octree.getBoundingBoxX1()[1];
      x1[2] = octree.getBoundingBoxX1()[2];
      x2[0] = octree.getBoundingBoxX2()[0];
      x2[1] = octree.getBoundingBoxX2()[1];
      x2[2] = octree.getBoundingBoxX2()[2];
      dx[0] = x2[0] - x1[0];
      dx[1] = x2[1] - x1[1];
      dx[2] = x2[2] - x1[2];
      //
      vector_t x;
      x[0] = x1[0] + m_Ix * dx[0] / scalar_t(65535);
      x[1] = x1[1] + m_Iy * dx[1] / scalar_t(65535);
      x[2] = x1[2] + m_Iz * dx[2] / scalar_t(65535);
      return x;
    };
    //
    std::uint8_t compare(const vertex_t& v) const
    {
      static const std::int8_t x = 0b00000001;
      static const std::int8_t y = 0b00000010;
      static const std::int8_t z = 0b00000100;
      //
      std::int8_t result = 0;
      if (m_Ix == v.m_Ix) result |= x;
      if (m_Iy == v.m_Iy) result |= y;
      if (m_Iz == v.m_Iz) result |= z;
      return result;
    }
    //
    bool operator==(const vertex_t& v) const
    {
      return (m_Ix == v.m_Ix && m_Iy == v.m_Iy && m_Iz == v.m_Iz && m_Index == v.m_Index);
    };
    //
    std::string toString() const
    {
      std::string s = "(";
      s += std::to_string(m_Ix) + ",";
      s += std::to_string(m_Iy) + ",";
      s += std::to_string(m_Iz) + ")";      
      return s;
    };
    //
    size_t hash() const
    {
      return (size_t(m_Ix) << 32) + (size_t(m_Iy) << 16) + size_t(m_Iz);
    };
    //
    bool operator<(const vertex_t& v) const
    {
      return hash() < v.hash();
    };
  };

  struct vertex_hash_t
  {
    size_t operator()(const vertex_t& v) const
    {
      return v.hash();
    }
  };

  friend class vertex_t;


  /**
    * @brief A node in the octree (i.e. a bucket or a parent node).
    * This is not to be confused with a geometric vertex in the octree, which is
    * described by the vertex_t type.
    */
  struct node_t
  {
    std::vector<int>      m_ItemIndices;
    vector_t              m_X1;
    vector_t              m_X2;
    vector_t              m_Centre;
    scalar_t              m_Diagonal;
    std::vector<vertex_t> m_Vertices;    
    int                   m_Level = 0;
    int                   m_Parent = -1;

    node_t() {};

    node_t(const Octree<VEC,ITEM,CHECK>& octree, vertex_t v1, vertex_t v2, int parent) : m_Parent(parent)
    {
      if (parent == -1) {
        m_Level = 0;
      } else {
        m_Level = octree.m_Nodes[parent].m_Level + 1;
      }
      vector_t x1 = v1.getCoordinates(octree);
      vector_t x2 = v2.getCoordinates(octree);
      vector_t dx = x2 - x1;
      //
      m_Centre   = 0.5 * (x1 + x2);
      m_X1       = m_Centre - 0.51*dx;
      m_X2       = m_Centre + 0.51*dx;
      m_Diagonal = (x2 - x1).abs();
      //
      // store or create the vertices!!
      //
      m_Vertices.resize(8);
      m_Vertices[0] = v1;
      m_Vertices[1] = vertex_t(v2.m_Ix, v1.m_Iy, v1.m_Iz);
      m_Vertices[2] = vertex_t(v1.m_Ix, v2.m_Iy, v1.m_Iz);
      m_Vertices[3] = vertex_t(v2.m_Ix, v2.m_Iy, v1.m_Iz);
      m_Vertices[4] = vertex_t(v1.m_Ix, v1.m_Iy, v2.m_Iz);
      m_Vertices[5] = vertex_t(v2.m_Ix, v1.m_Iy, v2.m_Iz);
      m_Vertices[6] = vertex_t(v1.m_Ix, v2.m_Iy, v2.m_Iz);
      m_Vertices[7] = v2;
    };

    std::vector<int> getItems(const Octree<VEC,ITEM,CHECK> &octree) const 
    {
      std::vector<int> items;
      if (m_Children.size() == 0) {
        if (m_ItemIndices.size() == 0) {
          auto &parent = octree.m_Nodes[m_Parent];
          for (int child_index : parent.m_Children) {
            const node_t &child = octree.m_Nodes[child_index];
            if (child.m_ItemIndices.size() != 0) {
              items.insert(items.end(), child.m_ItemIndices.begin(), child.m_ItemIndices.end());
            }
          }
        } else {
          items = m_ItemIndices;
        }
      }
      return items;
    };

    std::string toString() const
    {
      std::string s = "(";
      s += m_Vertices[0].toString() + " ";
      s += m_Vertices[7].toString();
      s += " level:" + std::to_string(m_Level);
      s += " children:" + std::to_string(m_Children.size());
      s += " items:" + std::to_string(m_ItemIndices.size());
      s += ")";
      return s;
    };

    std::vector<int> m_Children;
  };


protected: // attributes

  std::vector<node_t>           m_Nodes;
  std::vector<item_t>           m_Items;
  int                           m_MaxBucketSize;
  bool                          m_ApproximateSearch = false;
  bool                          m_AllowEmptyBuckets = false;
  vector_t                      m_X1;
  vector_t                      m_X2;
  bool                          m_CustomBoundingBox = false;
  std::vector<std::vector<int>> m_Vertex2Node;
  std::vector<std::vector<int>> m_SearchList;
  int                           m_MaxLevel = 0;  
  static const int              m_MaxLevelLimit = 15;


protected: // methods

  int find(vector_t x, int node_index) const
  {
    if (m_Nodes[node_index].m_Children.size() == 0) {
      scalar_t min_dist = 10*m_Nodes[node_index].m_Diagonal;
      int closest_item_index = -1;
      //
      // search all items in the search list
      //
      for (int item_index : m_SearchList[node_index]) {
        vector_t x1, x2;
        check_t::boundingBox(m_Items[item_index], x1, x2);
        vector_t x_item = 0.5 * (x1 + x2);
        scalar_t dist = (x - x_item).abs();
        if (dist < min_dist) {
          min_dist = dist;
          closest_item_index = item_index;
        }
      }
      return closest_item_index;
    }
    //
    if (x[0] > m_Nodes[node_index].m_Centre[0]) {
      if (x[1] > m_Nodes[node_index].m_Centre[1]) {
        if (x[2] > m_Nodes[node_index].m_Centre[2]) {
          return find(x, m_Nodes[node_index].m_Children[7]);
        } else {
          return find(x, m_Nodes[node_index].m_Children[3]);
        }
      } else {
        if (x[2] > m_Nodes[node_index].m_Centre[2]) {
          return find(x, m_Nodes[node_index].m_Children[5]);
        } else {
          return find(x, m_Nodes[node_index].m_Children[1]);
        }
      }
    } else {
      if (x[1] > m_Nodes[node_index].m_Centre[1]) {
        if (x[2] > m_Nodes[node_index].m_Centre[2]) {
          return find(x, m_Nodes[node_index].m_Children[6]);
        } else {
          return find(x, m_Nodes[node_index].m_Children[2]);
        }
      } else {
        if (x[2] > m_Nodes[node_index].m_Centre[2]) {
          return find(x, m_Nodes[node_index].m_Children[4]);
        } else {
          return find(x, m_Nodes[node_index].m_Children[0]);
        }
      }
    }
    return -1;
  }

  int findBucket(vector_t x, int node_index) const
  {
    if (node_index == 0) {
      if (!isInsideCartesianBox(x, m_X1, m_X2)) {
        return -1;
      }
    }
    if (m_Nodes[node_index].m_Children.size() == 0) {
      return node_index;
    }
    //
    if (x[0] > m_Nodes[node_index].m_Centre[0]) {
      if (x[1] > m_Nodes[node_index].m_Centre[1]) {
        if (x[2] > m_Nodes[node_index].m_Centre[2]) {
          return findBucket(x, m_Nodes[node_index].m_Children[7]);
        } else {
          return findBucket(x, m_Nodes[node_index].m_Children[3]);
        }
      } else {
        if (x[2] > m_Nodes[node_index].m_Centre[2]) {
          return findBucket(x, m_Nodes[node_index].m_Children[5]);
        } else {
          return findBucket(x, m_Nodes[node_index].m_Children[1]);
        }
      }
    } else {
      if (x[1] > m_Nodes[node_index].m_Centre[1]) {
        if (x[2] > m_Nodes[node_index].m_Centre[2]) {
          return findBucket(x, m_Nodes[node_index].m_Children[6]);
        } else {
          return findBucket(x, m_Nodes[node_index].m_Children[2]);
        }
      } else {
        if (x[2] > m_Nodes[node_index].m_Centre[2]) {
          return findBucket(x, m_Nodes[node_index].m_Children[4]);
        } else {
          return findBucket(x, m_Nodes[node_index].m_Children[0]);
        }
      }
    }
    return -1;
  }

  /**
    * create unique vertex list and set the index field in the vertices
    */
  void buildVertexStructures() 
  {
    std::unordered_map<uint64_t, uint32_t> vertex_map;
    for (int i = 0; i < m_Nodes.size(); ++i) {
      for (int j = 0; j < m_Nodes[i].m_Vertices.size(); ++j) {
        vertex_t& v = m_Nodes[i].m_Vertices[j];
        uint64_t key = (uint64_t(v.m_Ix) << 32) + (uint64_t(v.m_Iy) << 16) + uint64_t(v.m_Iz);
        auto it = vertex_map.find(key);
        if (it == vertex_map.end()) {
          v.m_Index = vertex_map.size();
          vertex_map[key] = v.m_Index;
        } else {
          v.m_Index = it->second;
        }
      }
    }
    m_Vertex2Node.resize(vertex_map.size());
    for (int i = 0; i < m_Nodes.size(); ++i) {
      for (int j = 0; j < m_Nodes[i].m_Vertices.size(); ++j) {
        vertex_t v = m_Nodes[i].m_Vertices[j];
        m_Vertex2Node[v.m_Index].push_back(i);
      }
    }
  }

  std::vector<vertex_t> getUniqueVertices(int node_index) const
  {
    std::unordered_set<vertex_t, vertex_hash_t> vertices;
    scalar_t radius = 1.5*m_Nodes[node_index].m_Diagonal;
    vector_t x_node = m_Nodes[node_index].m_Centre;
    for (auto vertex : m_Nodes[node_index].m_Vertices) {
      vertices.insert(vertex);
    }
    //
    int num_new_vertices;
    do {
      int num_vertices = vertices.size();
      std::vector<vertex_t> old_vertices(vertices.begin(), vertices.end());
      //
      for (auto old_vertex : old_vertices) {
        for (int node_index : m_Vertex2Node[old_vertex.m_Index]) {
          for (auto vertex : m_Nodes[node_index].m_Vertices) {
            auto x_vertex = vertex.getCoordinates(*this);
            if ((x_vertex - x_node).abs() <= radius) {
              vertices.insert(vertex);
            }
          }
        }
      }
      num_new_vertices = vertices.size() - num_vertices;
    } while (num_new_vertices > 0);
    //
    std::vector<vertex_t> result(vertices.begin(), vertices.end());
    return result;
  }

  std::vector<int> getSearchItems(vector_t x_search, std::vector<vertex_t> vertices, scalar_t radius)
  {
    std::unordered_set<int> items;
    //
    for (auto vertex : vertices) {
      for (int node_index : m_Vertex2Node[vertex.m_Index]) {
        if (m_Nodes[node_index].m_Children.size() == 0) {
          auto node_items = m_Nodes[node_index].getItems(*this);
          for (int item_index : node_items) {
            vector_t x1, x2;
            check_t::boundingBox(m_Items[item_index], x1, x2);
            if ((x1 - x_search).abs() <= radius) {
              items.insert(item_index);
            } else if ((x2 - x_search).abs() <= radius) {
              items.insert(item_index);
            }
          }
        }
      }
    }
    //
    std::vector<int> result(items.begin(), items.end());
    return result;
  }


public:

  Octree(int max_bucket_size=10) : m_MaxBucketSize(max_bucket_size) {};

  /**
    * @brief Get the first corner of the overall bounding box.
    * @return The first corner of the overall bounding box.
    */
  const vector_t& getBoundingBoxX1() const { return m_X1; };

  /**
    * @brief Get the second corner of the overall bounding box.
    * @return The second corner of the overall bounding box.
    */
  const vector_t& getBoundingBoxX2() const { return m_X2; };

  void setMaximalBucketSize(int max_bucket_size) { m_MaxBucketSize = max_bucket_size; }
  void approximateSearchOn() { m_ApproximateSearch = true; }
  void approximateSearchOff() { m_ApproximateSearch = false; }
  void allowEmptyBucketsOn() { m_AllowEmptyBuckets = true; }
  void allowEmptyBucketsOff() { m_AllowEmptyBuckets = false; }

  void setBoundingBox(vector_t x1, vector_t x2) 
  {
    m_X1 = x1;
    m_X2 = x2; 
    m_CustomBoundingBox = true;
  }

  template <typename C>
  void setItems(const C items)
  {
    m_Nodes.clear();
    if (!m_CustomBoundingBox) {
      bool first_item = true;
      for (auto item : items) {
        vector_t x1, x2;
        check_t::boundingBox(item, x1, x2);
        if (first_item) {
          m_X1 = x1;
          m_X2 = x2;
          first_item = false;
        } else {
          for (int i = 0; i < x1.size(); ++i) {
            m_X1[i] = std::min(m_X1[i], x1[i]);
            m_X2[i] = std::max(m_X2[i], x2[i]);
          }
        }
      }
      vector_t x_centre = 0.50 * (m_X1 + m_X2);
      vector_t Dx       = 0.55 * (m_X2 - m_X1);
      for (int i = 0; i < x_centre.size(); ++i) {
        Dx[i] = std::max(Dx[0], std::max(Dx[1], Dx[2]));
      }
      m_X1 = x_centre - Dx;
      m_X2 = x_centre + Dx;
    }
    vertex_t v1, v2;
    v2.m_Ix = 65535;
    v2.m_Iy = 65535;
    v2.m_Iz = 65535;
    node_t root(*this, v1, v2, -1);
    m_Nodes.push_back(root);
    //
    // copy all items into the root node and the coordinates into m_Items
    //
    m_Nodes[0].m_ItemIndices.resize(items.size());
    m_Items.resize(items.size());
    for (int i = 0; i < items.size(); ++i) {
      m_Nodes[0].m_ItemIndices[i] = i;
      m_Items[i] = items[i];
    }
    //
    // recursively divide the root node
    //
    m_MaxLevel = 0;
    refineNode(0);
    buildVertexStructures();
    //
    // create the search lists
    //
    m_SearchList.resize(m_Nodes.size());
    for (int i = 0; i < m_Nodes.size(); ++i) {
      if (m_Nodes[i].m_Children.size() == 0) {
        if (m_ApproximateSearch) {
          if (m_Nodes[i].m_ItemIndices.size() == 0 && !m_AllowEmptyBuckets) {
            auto& parent = m_Nodes[m_Nodes[i].m_Parent];
            for (int child_index : parent.m_Children) {
              const node_t& child = m_Nodes[child_index];
              if (child.m_ItemIndices.size() != 0) {
                m_SearchList[i].insert(m_SearchList[i].end(), child.m_ItemIndices.begin(), child.m_ItemIndices.end());
              }
            }
          } else {
            m_SearchList[i] = m_Nodes[i].m_ItemIndices;
          }
        } else {
          auto vertices = getUniqueVertices(i);
          m_SearchList[i] = getSearchItems(m_Nodes[i].m_Centre, vertices, 1.5*m_Nodes[i].m_Diagonal);
        }
      }
    }
  }

  void refineNode(int node_index)
  {
    //std::cout << "refineNode " << node_index << "\n  -" << m_Nodes.size() << " nodes in total\n  - level " << m_Nodes[node_index].m_Level << std::endl;
    vector_t x1 = m_Nodes[node_index].m_X1;
    vector_t x2 = m_Nodes[node_index].m_X2;
    vector_t xc = 0.5 * (x1 + x2);
    auto lx = x2[0] - x1[0];
    auto ly = x2[1] - x1[1];
    auto lz = x2[2] - x1[2];
    if (m_Nodes[node_index].m_ItemIndices.size() <= m_MaxBucketSize) {
      return;
    }
    if (m_Nodes[node_index].m_Level == m_MaxLevelLimit) {
      return;
    }
    if (m_Nodes[node_index].m_Level > m_MaxLevelLimit) {
      EDL_BUG;
    }
    //
    const std::vector<vertex_t>& vertices = m_Nodes[node_index].m_Vertices;
    const vertex_t& v1 = vertices[0];
    const vertex_t& v2 = vertices[7];
    vertex_t vc(v1, v2);
    if (v1.compare(vc) != 0 || v2.compare(vc) != 0) {
      std::cout << "ERROR: vertex comparison failed" << std::endl;
      std::cout << "v1 = " << v1.toString() << std::endl;
      std::cout << "v2 = " << v2.toString() << std::endl;
      EDL_BUG;
    }
    //
    // distribute the items to the potential children
    //
    std::vector<node_t> potential_children(8);
    uint16_t ix1 = v1.m_Ix;
    uint16_t iy1 = v1.m_Iy;
    uint16_t iz1 = v1.m_Iz;
    uint16_t ix2 = v2.m_Ix;
    uint16_t iy2 = v2.m_Iy;
    uint16_t iz2 = v2.m_Iz;
    potential_children[0] = node_t(*this, vertices[0], vc, node_index);
    potential_children[1] = node_t(*this, vertex_t(vertices[0], vertices[1]), vertex_t(vertices[1], vertices[7]), node_index);
    potential_children[2] = node_t(*this, vertex_t(vertices[0], vertices[2]), vertex_t(vertices[2], vertices[7]), node_index);
    potential_children[3] = node_t(*this, vertex_t(vertices[0], vertices[3]), vertex_t(vertices[3], vertices[7]), node_index);
    potential_children[4] = node_t(*this, vertex_t(vertices[0], vertices[4]), vertex_t(vertices[4], vertices[7]), node_index);
    potential_children[5] = node_t(*this, vertex_t(vertices[0], vertices[5]), vertex_t(vertices[5], vertices[7]), node_index);
    potential_children[6] = node_t(*this, vertex_t(vertices[0], vertices[6]), vertex_t(vertices[6], vertices[7]), node_index);
    potential_children[7] = node_t(*this, vc, vertices[7], node_index);
    //
    for (int i = 0; i < 8; ++i) {
      potential_children[i].m_ItemIndices.reserve(m_Nodes[node_index].m_ItemIndices.size()/8);
    }
    for (int item_index : m_Nodes[node_index].m_ItemIndices) {
      item_t p = m_Items[item_index];
      for (int i = 0; i < 8; ++i) {
        if (check_t::isInsideCartesianBox(p, potential_children[i].m_X1, potential_children[i].m_X2)) {
          potential_children[i].m_ItemIndices.push_back(item_index);
        }
      }
    }
    //
    // create the children
    //
    m_Nodes[node_index].m_Children.resize(8);
    for (int i = 0; i < 8; ++i) {
      m_Nodes.push_back(potential_children[i]);
      m_MaxLevel = std::max(m_MaxLevel, m_Nodes.back().m_Level);
      for (int j = 0; j < potential_children[i].m_ItemIndices.size(); ++j) {
        int item_index = potential_children[i].m_ItemIndices[j];
        if (item_index == -1) {
          EDL_BUG;
        }
      }
      m_Nodes[node_index].m_Children[i] = m_Nodes.size() - 1;
    }
    m_Nodes[node_index].m_ItemIndices.clear();
    //
    // recursively refine the children
    //
    for (int i_child = 0; i_child < m_Nodes[node_index].m_Children.size(); ++i_child) {
      int child_index = m_Nodes[node_index].m_Children[i_child];
      refineNode(child_index);
    }
  }

  int nearestItemIndex(vector_t x) const
  {
    return find(x, 0);
  }

  std::vector<int> getBucketContents(vector_t x) const
  {
    std::vector<int> items;
    int node_index = findBucket(x, 0);
    if (node_index != -1) {
      items = m_Nodes[node_index].m_ItemIndices;
    }
    return items;
  }

  int findMatchingItem(const vector_t& x, std::vector<int> items = {}) const
  {
    if (items.empty()) {
      items = getBucketContents(x);
    }
    for (int item_index : items) {
      if (check_t::match(x, m_Items[item_index])) {
        return item_index;
      }
    }
    return -1;
  }

  void dbgPrint()
  {
    for (int i = 0; i < m_Nodes.size(); ++i) {
      const node_t& node = m_Nodes[i];
      if (node.m_Children.size() == 0) {
        bool print = false;
        for (int idx : node.m_ItemIndices) {
          if (idx == 748) {
            print = true;
          }
        }
        if (print) {
          std::cout << "bucket " << i << " : ";
          std::cout << node.m_X1 << " ";
          std::cout << node.m_X2 << " : ";
          for (int idx : node.m_ItemIndices) {
            std::cout << idx << " ";
          }
          std::cout << std::endl;
        }
      }
    }
  }

  int depth() const
  {
    return m_MaxLevel;
  }

  /**
    * @brief Write the octree to a legacy VTK file.
    * This method does not require the VTK library.
    * @param file_name The name of the file to write.
    */
  void writeVtkFile(std::string file_name)
  {
    std::ofstream file(file_name);
    file << "# vtk DataFile Version 2.0" << std::endl;
    file << "Octree" << std::endl;
    file << "ASCII" << std::endl;
    file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    //
    // first count the number of vertices
    //
    buildVertexStructures();
    std::map<int,vector_t> vertex_map;
    for (auto node : m_Nodes) {
      if (node.m_Children.size() == 0) {
        for (auto v : node.m_Vertices) {
          vertex_map[v.m_Index] = v.getCoordinates(*this);
        }
      }
    }
    //
    // write the vertices
    //
    file << "POINTS " << vertex_map.size() << " float" << std::endl;
    for (auto V : vertex_map) {
      auto x = V.second;
      file << x[0] << " " << x[1] << " " << x[2] << std::endl;
    }
    //
    // write the cells
    //
    int num_cells = 0;
    for (auto node : m_Nodes) {
      if (node.m_Children.size() == 0) {
        num_cells += 1;
      }
    }
    file << "CELLS " << num_cells << " " << 9*num_cells << std::endl;
    for (int i = 0; i < m_Nodes.size(); ++i) {
      if (m_Nodes[i].m_Children.size() == 0) {
        file << "8 ";
        for (auto vertex : m_Nodes[i].m_Vertices) {
          file << vertex.m_Index << " ";
        }
        file << std::endl;
      }
    }
    //
    // write the cell types
    //
    file << "CELL_TYPES " << num_cells << std::endl;
    for (int i = 0; i < num_cells; ++i) {
      file << "11" << std::endl;
    }
    //
    // write the data
    //
    file << "CELL_DATA " << num_cells << std::endl;
    file << "SCALARS index int 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < m_Nodes.size(); ++i) {
      if (m_Nodes[i].m_Children.size() == 0) {
        file << i << std::endl;
      }
    }
    file << "SCALARS level int 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < m_Nodes.size(); ++i) {
      if (m_Nodes[i].m_Children.size() == 0) {
        file << m_Nodes[i].m_Level << std::endl;
      }
    }
    file << "SCALARS items int 1" << std::endl;
    file << "LOOKUP_TABLE default" << std::endl;
    for (int i = 0; i < m_Nodes.size(); ++i) {
      if (m_Nodes[i].m_Children.size() == 0) {
        file << m_Nodes[i].m_ItemIndices.size() << std::endl;
      }
    }
  }

};

} // namespace EDL_NAMESPACE


// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

#include <random>
#include "edl/geometrytools.h"

TEST_CASE("Octree__regularly_spaced_items")
{
  //
  // This doesn't test much...
  // Success is defined as it doesn't crash and finds an item.
  // If the item is correct or not is not checked.
  // It is mainly intended to check things during development
  //
  using namespace edl;
  using std::vector;
  typedef float real;
  typedef MathVector<StaticVector<real,3>> vec_t;
  typedef Octree<vec_t, vec_t, OctreeVectorCheck<vec_t>> octree_t;
  //
  int N = 5;
  vector<vec_t> items(N*N*N);
  real delta = 1.0 / (N-1);
  for (int i = 0; i < N; ++i) {
    real x = i * delta;
    for (int j = 0; j < N; ++j) {
      real y = j * delta;
      for (int k = 0; k < N; ++k) {
        real z = k * delta;
        items[i*N*N + j*N + k] = vec_t(x,y,z);
      }
    }
  }
  //
  // create the octree
  //
  octree_t octree(10);
  octree.setItems(items);
  //
  // create a set of test items
  //
  vector<vec_t> test_items(10);
  delta = 1.0 / (test_items.size() - 1);
  for (int i = 0; i < test_items.size(); ++i) {
    real x = i * delta;
    test_items[i] = vec_t(x,x,x);
  }
  //
  // perform the tests
  //
  for (int i = 0; i < test_items.size(); ++i) {
    vec_t p = test_items[i];
    int   j = octree.nearestItemIndex(p);
    CHECK(j >= 0);
  }
}

TEST_CASE("Octree__random_Items_search")
{
  using namespace edl;
  using std::vector; 
  typedef float real;
  typedef MathVector<StaticVector<real,3>> vec_t;
  typedef Octree<vec_t, vec_t, OctreeVectorCheck<vec_t>> octree_t;
  //
  const int  num_Items = 20000;
  const int  num_tests = 1000;
  const real range     = 3;
  QuickTimer timer;
  //
  // create a set of random items on a sphere with radius range/3
  //
  std::mt19937 gen(42);
  std::uniform_real_distribution<real> dist1(-M_PI, M_PI);
  vector<vec_t> items(num_Items);
  vec_t y_axis(0,1,0);
  vec_t z_axis(0,0,1);
  for (int i = 0; i < num_Items; ++i) {
    items[i] = vec_t(range/3,0,0);
    items[i] = rotate(items[i], y_axis, dist1(gen));
    items[i] = rotate(items[i], z_axis, dist1(gen));
  }
  //
  // create the octree
  //
  timer.start();
  octree_t octree(100);
  octree.setBoundingBox(vec_t(-range/2,-range/2,-range/2), vec_t(range/2,range/2,range/2));
  octree.setItems(items);  
  timer.stop();
  std::cout << "Octree construction time : " << timer.milliseconds() << " ms" << std::endl;
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
  for (int i = 0; i < test_items.size(); ++i) {
    vec_t p = test_items[i];
    real  d_min = 1e10;
    for (int j = 0; j < num_Items; ++j) {
      vec_t q = items[j];
      real  d = (p - q).abs();
      if (d < d_min) {
        d_min = d;
        nearest_item[i] = j;
      }
    }
  }
  timer.stop();
  std::cout << "Brute force search time  : " << timer.milliseconds() << " ms" << std::endl;
  //
  // perform the tests
  //
  timer.restart();
  for (int i = 0; i < test_items.size(); ++i) {
    vec_t p = test_items[i];
    int   j = octree.nearestItemIndex(p);
    CHECK(j == nearest_item[i]);
  }
  timer.stop();
  std::cout << "Octree search time       : " << timer.milliseconds() << " ms" << std::endl;
}

TEST_CASE("Octree__random_Items_search_approximate")
{
  using namespace edl;
  using std::vector;
  typedef float real;
  typedef MathVector<StaticVector<real,3>> vec_t;
  typedef Octree<vec_t, vec_t, OctreeVectorCheck<vec_t>> octree_t;
  //
  const int  num_Items  = 2000;
  const int  num_tests   = 100;
  const real range       = 3;
  //
  // create a set of random items
  //
  //random_device rd;
  std::mt19937 gen(42);
  std::uniform_real_distribution<real> dist1(-range/2, range/2);
  vector<vec_t> items(num_Items);
  for (int i = 0; i < num_Items; ++i) {
    items[i] = vec_t(dist1(gen), dist1(gen), dist1(gen));
  }
  //
  // create the octree
  //
  octree_t octree(10);
  octree.approximateSearchOn();
  octree.setItems(items);
  //
  // create a set of test items
  //
  vector<vec_t> test_items(num_tests);
  for (int i = 0; i < num_tests; ++i) {
    test_items[i] = vec_t(dist1(gen), dist1(gen), dist1(gen));
  }
  //
  // perform the tests
  //
  for (int i = 0; i < num_tests; ++i) {
    vec_t p = test_items[i];
    int   j = octree.nearestItemIndex(p);
    //
    // check if we find any item at all
    // We cannot check if the item is correct, because the search is approximate
    //
    CHECK(j >= 0);
  }
}

// test a search for tetrahedra

namespace OctreeTetraSearch 
{

  using namespace edl;
  using std::vector;
  typedef float real;
  typedef MathVector<StaticVector<real,3>> vec_t;

  struct tetra_t
  {
    vector<vec_t> x;
    tetra_t() : x(4) {};
  };

  struct TetraCheck
  {    
    static bool isInsideCartesianBox(const tetra_t& tetra, const vec_t& x1, const vec_t& x2)
    {
      std::vector<vec_t> tetra_vertices(4);
      tetra_vertices[0] = tetra.x[0];
      tetra_vertices[1] = tetra.x[1];
      tetra_vertices[2] = tetra.x[2];
      tetra_vertices[3] = tetra.x[3];
      return tetraIsInsideCartesianBox(tetra_vertices, x1, x2);
    }

    static bool match(const vec_t& x, const tetra_t& tetra)
    {
      return isPointInTetra(tetra.x[0], tetra.x[1], tetra.x[2], tetra.x[3], x);
    }

    static void boundingBox(const tetra_t& tetra, vec_t& x1, vec_t& x2)
    {
      x1 = tetra.x[0];
      x2 = tetra.x[0];
      for (int i = 1; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
          x1[j] = std::min(x1[j], tetra.x[i][j]);
          x2[j] = std::max(x2[j], tetra.x[i][j]);
        }
      }
    }
  };

  typedef Octree<vec_t, tetra_t, TetraCheck> octree_t;

} // namespace OctreeTetraSearch

TEST_CASE("Octree_tetra_box_check")
{
  using namespace edl;
  using namespace edl;
  using namespace OctreeTetraSearch;
  //
  tetra_t tetra;
  tetra.x[0] = vec_t(0,0,0);
  tetra.x[1] = vec_t(1,0,0);
  tetra.x[2] = vec_t(0,1,0);
  tetra.x[3] = vec_t(0,0,1);
  //
  vec_t x1, x2;
  //
  x1 = vec_t( 2, 2, 2);
  x2 = vec_t( 3, 3, 3);
  //
  CHECK(!TetraCheck::isInsideCartesianBox(tetra, x1, x2));
}

TEST_CASE("Octree_tetra_search")
{
  using std::vector; using std::uint32_t;
  using namespace edl;
  using namespace OctreeTetraSearch;
  using OctreeTetraSearch::real;
  //
  std::string file_name = "../../tests/sphere_test.mesh";
  //string file_name = "../../tests/test_box.mesh";
  vector<vec_t>       points;
  vector<tetra_t>     tetras;
  vector<vector<int>> tetra_vertices;
  //
  // read points and tetras from file
  //
  std::cout << file_name << std::endl;
  std::ifstream file(file_name);
  CHECK(file.is_open());
  //
  int N;
  file >> N;
  points.resize(N);
  for (int i = 0; i < N; ++i) {
    file >> points[i][0] >> points[i][1] >> points[i][2];
  }
  //
  file >> N;
  tetras.resize(N);
  tetra_vertices.resize(N);
  uint32_t i_min = std::numeric_limits<uint32_t>::max();
  uint32_t i_max = 0;
  for (int i = 0; i < N; ++i) {
    tetra_vertices[i].resize(4);
    int dummy;
    file >> dummy;
    for (int j = 0; j < 4; ++j) {
      uint32_t idx;
      file >> idx;
      idx -= 1;
      tetras[i].x[j] = points[idx];
      tetra_vertices[i][j] = idx;
      i_min = std::min(i_min, idx);
      i_max = std::max(i_max, idx);
    }
  }
  CHECK(points.size() == 9221);
  CHECK(tetras.size() == 41604);
  //
  // compute the maximal number of tetras a vertex is part of
  //
  vector<int> vertex_count(points.size(), 0);
  for (auto tetra : tetra_vertices) {
    for (auto idx : tetra) {
      vertex_count[idx] += 1;
    }
  }
  int max_vertex_count = 0;
  for (auto count : vertex_count) {
    max_vertex_count = std::max(max_vertex_count, count);
  }
  //
  // write the tetra test mesh to a VTK file
  //
  /*
  ofstream vtk_file("octree_tetra_test_tetras.vtk");
  vtk_file << "# vtk DataFile Version 2.0" << endl;
  vtk_file << "Tetras" << endl;
  vtk_file << "ASCII" << endl;
  vtk_file << "DATASET UNSTRUCTURED_GRID" << endl;
  vtk_file << "POINTS " << points.size() << " float" << endl;
  for (auto p : points) {
    vtk_file << p[0] << " " << p[1] << " " << p[2] << endl;
  }
  vtk_file << "CELLS " << tetras.size() << " " << 5*tetras.size() << endl;
  for (auto tetra : tetra_vertices) {
    vtk_file << "4 ";
    for (auto i : tetra) {
      vtk_file << i << " ";
    }
    vtk_file << endl;
  }
  vtk_file << "CELL_TYPES " << tetras.size() << endl;
  for (int i = 0; i < tetras.size(); ++i) {
    vtk_file << "10" << endl;
  }
  // write the tetra index (item_index)
  vtk_file << "CELL_DATA " << tetras.size() << endl;
  vtk_file << "SCALARS tetra_index int 1" << endl;
  vtk_file << "LOOKUP_TABLE default" << endl;
  for (int i = 0; i < tetras.size(); ++i) {
    vtk_file << i << endl;
  }
  */
  //
  // build the octree
  //
  octree_t octree(max_vertex_count);
  octree.approximateSearchOn();
  octree.allowEmptyBucketsOn();
  octree.setItems(tetras);
  //octree.writeVtkFile("octree_tetra_test_octree.vtk");
  //
  // loop over all Cartesian coordinates and check if we find a tetra
  // for positions which are inside the mesh (R1 < r < R2)
  //
  real  R1 = 0.255;
  real  R2 = 0.495;
  int   n  = 10;
  vec_t x1 = vec_t(-R2, -R2, -R2);
  vec_t x2 = vec_t( R2,  R2,  R2);
  for (int i = 0; i < n; ++i) {
    real x = x1[0] + (x2[0] - x1[0])*i/(n-1);
    for (int j = 0; j < n; ++j) {
      real y = x1[1] + (x2[1] - x1[1])*j/(n-1);
      for (int k = 0; k < n; ++k) {
        real z = x1[2] + (x2[2] - x1[2])*k/(n-1);
        vec_t p(x,y,z);
        real r = p.abs();
        if (R1 < r && r < R2) {
          int idx = octree.findMatchingItem(p);
          CHECK(idx >= 0);
          CHECK(idx < tetras.size());
          CHECK(TetraCheck::match(p, tetras[idx]));
        }
      }
    }
  }
}

#endif // OCTREE_H