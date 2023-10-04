// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2021 enGits GmbH                                         +
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

#ifndef OCTREE_H
#define OCTREE_H

#include <bits/stdint-uintn.h>
#include <cstdint>
#include <iostream>

#include <ostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>

#include "edl/edl.h"
#include "edl/edlerror.h"
#include "edl/geometrytools.h"

namespace EDL_NAMESPACE
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

template <typename VEC>
class Octree
{

public: // data types

  typedef VEC vector_t;
  typedef typename vector_t::value_type scalar_t;


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
    vertex_t(vector_t x, const Octree<VEC>& octree)
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
    vector_t getCoordinates(const Octree<VEC>& octree) const
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
    uint64_t hash() const
    {
      return (uint64_t(m_Ix) << 32) + (uint64_t(m_Iy) << 16) + uint64_t(m_Iz);
    };
    //
    bool operator<(const vertex_t& v) const
    {
      return hash() < v.hash();
    };
  };

  friend class vertex_t;


  /**
    * @brief A node in the octree (i.e. a bucket or a parent node).
    * This is not to be confused with a geometric vertex in the octree, which is
    * described by the vertex_t type.
    */
  struct node_t
  {  
    std::vector<int>      m_PointIndices;
    vector_t              m_X1;
    vector_t              m_X2;
    vector_t              m_Centre;
    scalar_t              m_Diagonal;
    std::vector<vertex_t> m_Vertices;    
    int                   m_Level = 0;
    int                   m_Parent = -1;

    node_t() {};

    node_t(const Octree<VEC>& octree, vertex_t v1, vertex_t v2, int parent) : m_Parent(parent)
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

    std::vector<int> getPoints(const Octree<VEC> &octree) const 
    {
      std::vector<int> points;
      if (m_Children.size() == 0) {
        if (m_PointIndices.size() == 0) {
          auto &parent = octree.m_Nodes[m_Parent];
          for (int child_index : parent.m_Children) {
            const node_t &child = octree.m_Nodes[child_index];
            if (child.m_PointIndices.size() != 0) {
              points.insert(points.end(), child.m_PointIndices.begin(), child.m_PointIndices.end());
            }
          }
        } else {
          points = m_PointIndices;
        }
      }
      return points;
    };

    std::string toString() const
    {
      std::string s = "(";
      s += m_Vertices[0].toString() + " ";
      s += m_Vertices[7].toString();
      s += " level:" + std::to_string(m_Level);
      s += " children:" + std::to_string(m_Children.size());
      s += " points:" + std::to_string(m_PointIndices.size());
      s += ")";
      return s;
    };

    std::vector<int> m_Children;
  };


protected: // attributes

  std::vector<node_t>           m_Nodes;
  std::vector<vector_t>         m_Points;
  std::vector<vertex_t>         m_Vertices;
  int                           m_MaxBucketSize;
  bool                          m_ApproximateSearch = false;
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
      int closest_point_index = -1;
      //
      // search all points in the search list
      //
      for (int point_index : m_SearchList[node_index]) {
        scalar_t dist = (x - m_Points[point_index]).abs();
        if (dist < min_dist) {
          min_dist = dist;
          closest_point_index = point_index;
        }
      }
      return closest_point_index;
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
    std::set<vertex_t> vertices;
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

  std::vector<int> getSearchPoints(vector_t x_search, std::vector<vertex_t> vertices, scalar_t radius)
  {
    std::unordered_set<int> points;
    //
    for (auto vertex : vertices) {
      for (int node_index : m_Vertex2Node[vertex.m_Index]) {
        if (m_Nodes[node_index].m_Children.size() == 0) {
          auto node_points = m_Nodes[node_index].getPoints(*this);
          for (int point_index : node_points) {
            if ((m_Points[point_index] - x_search).abs() <= radius) {
              points.insert(point_index);
            }
          }
        }
      }
    }
    //
    std::vector<int> result(points.begin(), points.end());
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

  void setBoundingBox(vector_t x1, vector_t x2) 
  {
    m_X1 = x1;
    m_X2 = x2; 
    m_CustomBoundingBox = true;
  }

  template <typename C>
  void setPoints(const C points)
  {
    m_Nodes.clear();
    if (!m_CustomBoundingBox) {
      findBoundingBox(points, m_X1, m_X2);
    }
    vertex_t v1, v2;
    v2.m_Ix = 65535;
    v2.m_Iy = 65535;
    v2.m_Iz = 65535;
    node_t root(*this, v1, v2, -1);
    m_Nodes.push_back(root);
    //
    // copy all points into the root node and the coordinates into m_Points
    //
    m_Nodes[0].m_PointIndices.resize(points.size());
    m_Points.resize(points.size());
    for (int i = 0; i < points.size(); ++i) {
      m_Nodes[0].m_PointIndices[i] = i;
      m_Points[i] = points[i];
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
          if (m_Nodes[i].m_PointIndices.size() == 0) {
            auto& parent = m_Nodes[m_Nodes[i].m_Parent];
            for (int child_index : parent.m_Children) {
              const node_t& child = m_Nodes[child_index];
              if (child.m_PointIndices.size() != 0) {
                m_SearchList[i].insert(m_SearchList[i].end(), child.m_PointIndices.begin(), child.m_PointIndices.end());
              }
            }
          } else {
            m_SearchList[i] = m_Nodes[i].m_PointIndices;
          }
        } else {
          auto vertices = getUniqueVertices(i);
          m_SearchList[i] = getSearchPoints(m_Nodes[i].m_Centre, vertices, 1.5*m_Nodes[i].m_Diagonal);
        }
      }
    }
  }

  void refineNode(int node_index)
  {
    if (m_Nodes[node_index].m_PointIndices.size() <= m_MaxBucketSize) {
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
    // distribute the points to the potential children
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
    for (int point_index : m_Nodes[node_index].m_PointIndices) {
      vector_t p = m_Points[point_index];
      for (int i = 0; i < 8; ++i) {
        if (isInsideCartesianBox(p,potential_children[i].m_X1, potential_children[i].m_X2, 0)) {
          // if (point_index == 314) {
          //   std::cout << "point 314 (" << p << ") is put in bucket ";
          //   std::cout << potential_children[i].toString() << std::endl;
          // }
          potential_children[i].m_PointIndices.push_back(point_index);
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
      for (int j = 0; j < potential_children[i].m_PointIndices.size(); ++j) {
        int point_index = potential_children[i].m_PointIndices[j];
        if (point_index == -1) {
          EDL_BUG;
        }
      }
      m_Nodes[node_index].m_Children[i] = m_Nodes.size() - 1;
    }
    m_Nodes[node_index].m_PointIndices.clear();
    //
    // recursively refine the children
    //
    for (int i_child = 0; i_child < m_Nodes[node_index].m_Children.size(); ++i_child) {
      int child_index = m_Nodes[node_index].m_Children[i_child];
      refineNode(child_index);
    }
  }

  int nearestPointIndex(vector_t x) const
  {
    return find(x, 0);
  }

  void dbgPrint()
  {
    for (int i = 0; i < m_Nodes.size(); ++i) {
      const node_t& node = m_Nodes[i];
      if (node.m_Children.size() == 0) {
        bool print = false;
        for (int idx : node.m_PointIndices) {
          if (idx == 748) {
            print = true;
          }
        }
        if (print) {
          std::cout << "bucket " << i << " : ";
          std::cout << node.m_X1 << " ";
          std::cout << node.m_X2 << " : ";
          for (int idx : node.m_PointIndices) {
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

};

} // namespace EDL_NAMESPACE


// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

#include <random>


TEST_CASE("Octree__regularly_spaced_points")
{
  //
  // This doesn't test much...
  // Success is defined as it doesn't crash and finds a point.
  // If the point is correct or not is not checked.
  // It is mainly intended to check things during development
  //
  using namespace edl;
  using namespace std;
  typedef float real;
  typedef MathVector<StaticVector<real,3>> vec_t;
  //
  int N = 5;
  vector<vec_t> points(N*N*N);
  real delta = 1.0 / (N-1);
  for (int i = 0; i < N; ++i) {
    real x = i * delta;
    for (int j = 0; j < N; ++j) {
      real y = j * delta;
      for (int k = 0; k < N; ++k) {
        real z = k * delta;
        points[i*N*N + j*N + k] = vec_t(x,y,z);
      }
    }
  }
  //
  // create the octree
  //
  Octree<vec_t> octree(10);
  octree.setPoints(points);
  //
  // create a set of test points
  //
  vector<vec_t> test_points(10);
  delta = 1.0 / (test_points.size() - 1);
  for (int i = 0; i < test_points.size(); ++i) {
    real x = i * delta;
    test_points[i] = vec_t(x,x,x);
  }
  //
  // perform the tests
  //
  for (int i = 0; i < test_points.size(); ++i) {
    vec_t p = test_points[i];
    int   j = octree.nearestPointIndex(p);
    CHECK(j >= 0);
  }
}

TEST_CASE("Octree__random_points_search")
{
  using namespace edl;
  using namespace std;
  typedef float real;
  typedef MathVector<StaticVector<real,3>> vec_t;
  //
  const int  num_points  = 1000;
  const int  num_tests   = 100;
  const real range       = 3;
  //
  // create a set of random points
  //
  //random_device rd;
  mt19937 gen(42);
  uniform_real_distribution<real> dist1(-range/2, range/2);
  vector<vec_t> points(num_points);
  for (int i = 0; i < num_points; ++i) {
    points[i] = vec_t(dist1(gen), dist1(gen), dist1(gen));
  }
  //
  // create the octree
  //
  Octree<vec_t> octree(10);
  octree.setPoints(points);
  //
  // create a set of test points
  //
  vector<vec_t> test_points(num_tests);
  for (int i = 0; i < num_tests; ++i) {
    test_points[i] = vec_t(dist1(gen), dist1(gen), dist1(gen));
  }
  //
  // brute force loop to find the reference solution
  //
  vector<int> nearest_point(num_tests);
  for (int i = 0; i < num_tests; ++i) {
    vec_t p = test_points[i];
    real  d_min = 1e10;
    for (int j = 0; j < num_points; ++j) {
      vec_t q = points[j];
      real  d = (p - q).abs();
      if (d < d_min) {
        d_min = d;
        nearest_point[i] = j;
      }
    }
  }
  //
  // perform the tests
  //
  //octree.dbgPrint();
  for (int i = 0; i < num_tests; ++i) {
    vec_t p = test_points[i];
    int   j = octree.nearestPointIndex(p);
    CHECK(j == nearest_point[i]);
  }
}

TEST_CASE("Octree__random_points_search_approximate")
{
  using namespace edl;
  using namespace std;
  typedef float real;
  typedef MathVector<StaticVector<real,3>> vec_t;
  //
  const int  num_points  = 1000;
  const int  num_tests   = 100;
  const real range       = 3;
  //
  // create a set of random points
  //
  //random_device rd;
  mt19937 gen(42);
  uniform_real_distribution<real> dist1(-range/2, range/2);
  vector<vec_t> points(num_points);
  for (int i = 0; i < num_points; ++i) {
    points[i] = vec_t(dist1(gen), dist1(gen), dist1(gen));
  }
  //
  // create the octree
  //
  Octree<vec_t> octree(10);
  octree.approximateSearchOn();
  octree.setPoints(points);
  //
  // create a set of test points
  //
  vector<vec_t> test_points(num_tests);
  for (int i = 0; i < num_tests; ++i) {
    test_points[i] = vec_t(dist1(gen), dist1(gen), dist1(gen));
  }
  //
  // perform the tests
  //
  for (int i = 0; i < num_tests; ++i) {
    vec_t p = test_points[i];
    int   j = octree.nearestPointIndex(p);
    //
    // check if we find any point at all
    // We cannot check if the point is correct, because the search is approximate
    //
    CHECK(j >= 0);
  }
}


#endif