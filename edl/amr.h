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

#ifndef AMRINDEX_H
#define AMRINDEX_H

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "edl/mathvector.h"

namespace edl
{
  
template <typename TIndex>
class AMRIndex 
{
public: // static attributes

  static const uint8_t m_MaxLevel = sizeof(TIndex) * 8 - 1;
  static const TIndex  negative1  = TIndex(-1);


private: // data types

  typedef TIndex index_type;


private: // attributes

  index_type m_I;
  index_type m_J;
  index_type m_K;  
  index_type m_Level;


public:

  AMRIndex() : m_I(-1), m_J(0), m_K(0), m_Level(0) 
  {
  }

  AMRIndex(index_type i, index_type j, index_type k, index_type level=0) : m_I(i), m_J(j), m_K(k), m_Level(level) 
  {    
  }

  void reduce()
  {
    if (m_Level == 0) {
      return;
    }
    while ((m_I % 2 == 0) && (m_J % 2 == 0) && (m_K % 2 == 0) && (m_Level > 0)) {
      m_I /= 2;
      m_J /= 2;
      m_K /= 2;
      --m_Level;
    }
  }

  bool valid() const
  {
    return (m_I != negative1) && (m_J != negative1) && (m_K != negative1) && (m_Level != negative1);
  }

  bool outOfBounds() const
  {
    return (m_I == negative1) && (m_J == 0) && (m_K == 0) && (m_Level == 0);
  }

  void setOutOfBounds()
  {
    m_I = -1;
    m_J = 0;
    m_K = 0;
    m_Level = 0;
  }

  AMRIndex reduced() const
  {
    AMRIndex<index_type> node(*this);
    node.reduce();
    return node;
  }

  AMRIndex indexOnLevel(index_type level) const
  {
    AMRIndex<index_type> node(*this);
    while (node.level() < level) {
      node.incrLevel();
    }
    while (node.level() > level && node.level() > 0) {
      node.reduce();
    }
    if (node.level() != level) {
      throw std::runtime_error("AMRIndex: indexOnLevel: invalid level");
    }
    return node;
  }

  void incrLevel()
  {
    m_I *= 2;
    m_J *= 2;
    m_K *= 2;
    ++m_Level;
  }

  AMRIndex incrementedLevel() const
  {
    AMRIndex<index_type> node(*this);
    node.incrLevel();
    return node;
  }

  void incrI()
  {
    m_I += 1;
  }

  AMRIndex incrementedI() const
  {
    AMRIndex<index_type> node(*this);
    node.incrI();
    return node;
  }

  void incrJ()
  {
    m_J += 1;
  }

  AMRIndex incrementedJ() const
  {
    AMRIndex<index_type> node(*this);
    node.incrJ();
    return node;
  }

  void incrK()
  {
    m_K += 1;
  }

  AMRIndex incrementedK() const
  {
    AMRIndex<index_type> node(*this);
    node.incrK();
    return node;
  }

  void decrI()
  {
    m_I -= 1;
  }

  AMRIndex decrementedI() const
  {
    AMRIndex<index_type> node(*this);
    node.decrI();
    return node;
  }

  void decrJ()
  {
    m_J -= 1;
  }

  AMRIndex decrementedJ() const
  {
    AMRIndex<index_type> node(*this);
    node.decrJ();
    return node;
  }

  void decrK()
  {
    m_K -= 1;
  }

  AMRIndex decrementedK() const
  {
    AMRIndex<index_type> node(*this);
    node.decrK();
    return node;
  }

  bool operator==(const AMRIndex& other) const
  {
    //return (m_I == other.i()) && (m_J == other.j()) && (m_K == other.k()) && (m_Level == other.level());
    if (m_Level != other.level()) {
      return false;
    }
    if (m_I != other.m_I) {
      return false;
    }
    if (m_J != other.m_J) {
      return false;
    }
    if (m_K != other.m_K) {
      return false;
    }
    return true;
  }

  void operator+=(const AMRIndex& other)
  {
    index_type i2     = other.i();
    index_type j2     = other.j();
    index_type k2     = other.k();
    index_type level2 = other.level();
    while (m_Level < level2) {
      m_I *= 2;
      m_J *= 2;
      m_K *= 2;
      ++m_Level;
    }
    while (m_Level > level2) {
      i2 *= 2;
      j2 *= 2;
      k2 *= 2;
      ++level2;
    }
    m_I += i2;
    m_J += j2;
    m_K += k2;
  }

  AMRIndex operator+(const AMRIndex& other) const
  {
    AMRIndex<index_type> node(*this);
    node += other;
    return node;
  }

  void operator-=(const AMRIndex& other)
  {
    index_type i2     = other.i();
    index_type j2     = other.j();
    index_type k2     = other.k();
    index_type level2 = other.level();
    while (m_Level < level2) {
      m_I *= 2;
      m_J *= 2;
      m_K *= 2;
      ++m_Level;
    }
    while (m_Level > level2) {
      i2 *= 2;
      j2 *= 2;
      k2 *= 2;
      ++level2;
    }
    m_I -= i2;
    m_J -= j2;
    m_K -= k2;
  }

  AMRIndex operator-(const AMRIndex& other) const
  {
    AMRIndex<index_type> node(*this);
    node -= other;
    return node;
  }

  bool operator!=(const AMRIndex& other) const
  {
    return !(*this == other);
  }

  void operator*=(index_type factor)
  {
    m_I *= factor;
    m_J *= factor;
    m_K *= factor;
  }

  AMRIndex operator*(index_type factor) const
  {
    AMRIndex<index_type> node(*this);
    node *= factor;
    return node;
  }

  void operator/=(index_type factor)
  {
    while ((m_I % factor != 0) || (m_J % factor != 0) || (m_K % factor != 0)) {
      if (m_Level > m_MaxLevel) {
        throw std::runtime_error("AMRIndex: division by factor not possible");
      }
      m_I *= 2;
      m_J *= 2;
      m_K *= 2;
      ++m_Level;
    }
    m_I /= factor;
    m_J /= factor;
    m_K /= factor;
  }

  AMRIndex operator/(index_type factor) const
  {
    AMRIndex<index_type> node(*this);
    node /= factor;
    return node;
  }

  AMRIndex parent() const
  {
    if (m_Level == 0) {
      return invalid();
    }
    AMRIndex<index_type> node(*this);
    node.m_I /= 2;
    node.m_J /= 2;
    node.m_K /= 2;
    --node.m_Level;
    return node;
  }

  index_type i()     const { return m_I; }
  index_type j()     const { return m_J; }
  index_type k()     const { return m_K; }
  index_type level() const { return m_Level; }


public: // static methods 

  static AMRIndex<index_type> invalid()
  {
    return AMRIndex<index_type>(-1, -1, -1, -1);
  }

  static AMRIndex<index_type> out_of_bounds()
  {
    return AMRIndex<index_type>(-1, 0, 0, 0);
  }

};

template <typename TIndex>
std::ostream& operator<<(std::ostream& os, const AMRIndex<TIndex>& node)
{
  os << "AMRIndex(" << node.i() << ", " << node.j() << ", " << node.k() << ", " << node.level() << ")";
  return os;
}

template <typename TIndex>
AMRIndex<TIndex> operator*(TIndex factor, const AMRIndex<TIndex>& node)
{
  return node * factor;
}


// ============================================================================

template <typename TIndex, typename Tijk, typename TVector>
class AMRMesh 
{
public: // static attributes



public: // data types

  typedef TIndex             index_type;
  typedef TVector            vector_type;
  typedef Tijk               ijk_type;
  typedef AMRIndex<ijk_type> amr_index_type;

  struct cell_t
  {
    index_type     linear_index;
    amr_index_type first_child;
    bool           is_outside = false;
  };

  struct node_t
  {
    index_type linear_index;
  };

  struct face_t
  {
    std::vector<amr_index_type> nodes;
    //
    index_type     cell1;
    index_type     cell2;
    bool           is_boundary = false;
    //
    face_t() {}
    face_t(const face_t& other) : nodes(other.nodes), cell1(other.cell1), cell2(other.cell2), is_boundary(other.is_boundary) {}
    face_t& operator=(const face_t& other)
    {
      nodes = other.nodes;
      cell1 = other.cell1;
      cell2 = other.cell2;
      is_boundary = other.is_boundary;
      return *this;
    }
    //
    bool operator<(const face_t& other) const
    {
      if (!is_boundary && other.is_boundary) {
        return true;
      }
      if (is_boundary && !other.is_boundary) {
        return false;
      }
      if (cell1 < other.cell1) {
        return true;
      }
      if (cell1 > other.cell1) {
        return false;
      }
      if (cell2 < other.cell2) {
        return true;
      }
      if (cell2 > other.cell2) {
        return false;
      }
      return false;
    }
    //
    void reverse()
    {
        // Reverse the nodes vector
        std::reverse(nodes.begin(), nodes.end());
        //
        // Swap cell1 and cell2
        std::swap(cell1, cell2);
    }
  };


private:; // attributes

  ijk_type    m_SizeI;
  ijk_type    m_SizeJ;
  ijk_type    m_SizeK;
  ijk_type    m_MaxLevel = 0;
  vector_type m_X1;
  vector_type m_X2;

  std::unordered_map<amr_index_type, cell_t> m_Cells;
  std::unordered_map<amr_index_type, node_t> m_Nodes;
  std::vector<face_t>                        m_Faces;

  
private: // methods

  void addNode(const amr_index_type& idx)
  {
    if (m_Nodes.find(idx.reduced()) == m_Nodes.end()) {
      node_t node;
      node.linear_index = m_Nodes.size();
      m_Nodes[idx.reduced()] = node;
    }
  }

  void addCell(const amr_index_type& idx)
  {
    m_MaxLevel = std::max(m_MaxLevel, idx.level());
    cell_t cell;
    cell.linear_index = m_Cells.size();
    cell.first_child  = amr_index_type();
    m_Cells[idx] = cell;
    auto node0 = idx;
    auto node1 = node0.incrementedI();
    auto node2 = node1.incrementedJ();
    auto node3 = node0.incrementedJ();
    auto node4 = node0.incrementedK();
    auto node5 = node1.incrementedK();
    auto node6 = node2.incrementedK();
    auto node7 = node3.incrementedK();
    addNode(node0);
    addNode(node1);
    addNode(node2);
    addNode(node3);
    addNode(node4);
    addNode(node5);
    addNode(node6);
    addNode(node7);
  }

  index_type sizeIOnLevel(index_type level) const
  {
    index_type size_i = m_SizeI;
    for (int i = 0; i < level; ++i) {
      size_i *= 2;
    }
    return size_i;
  }

  index_type sizeJOnLevel(index_type level) const
  {
    index_type size_j = m_SizeJ;
    for (int i = 0; i < level; ++i) {
      size_j *= 2;
    }
    return size_j;
  }

  index_type sizeKOnLevel(index_type level) const
  {
    index_type size_k = m_SizeK;
    for (int i = 0; i < level; ++i) {
      size_k *= 2;
    }
    return size_k;
  }

  amr_index_type cellNeighbour(const amr_index_type& idx1, const amr_index_type& idx2) const
  {
    auto idx2_it = m_Cells.find(idx2);
    if (idx2_it == m_Cells.end()) {
      if (idx2.level() > 0) {
        return cellNeighbour(idx1, idx2.parent());
      } else {
        return amr_index_type::out_of_bounds();
      }
    }
    cell_t cell2 = idx2_it->second;
    if (idx1.level() >= idx2.level() && !cell2.first_child.valid()) {
      return idx2;
    }
    return amr_index_type::invalid();
  }

  void completeFace(face_t& face) const
  {
    using namespace std;
    //
    vector<amr_index_type> all_nodes;
    all_nodes.reserve(8);
    for (int i = 0; i < face.nodes.size(); ++i) {
      auto idx0 = face.nodes[i];
      auto idx1 = face.nodes[(i + 1) % face.nodes.size()];
      all_nodes.push_back(idx0);
      auto nodes = nodesInBetween(idx0, idx1);
      for (auto idx : nodes) {
        all_nodes.push_back(idx);
      }
    }
    //
    face.nodes = all_nodes;
  }

  size_t numPotentialFaces()
  {
    size_t num_faces = 0;
    amr_index_type neighbours[6];
    auto leaf_cells = getLeafCellIndices();
    for (auto idx : leaf_cells) {
      auto cell = m_Cells.at(idx);
      neighbours[0] = cellNeighbourIM(idx);
      neighbours[1] = cellNeighbourIP(idx);
      neighbours[2] = cellNeighbourJM(idx);
      neighbours[3] = cellNeighbourJP(idx);
      neighbours[4] = cellNeighbourKM(idx);
      neighbours[5] = cellNeighbourKP(idx);
      //
      for (int i = 0; i < 6; ++i) {
        auto idx2 = neighbours[i];
        if (idx2.valid() || idx2.outOfBounds()) {
          bool use_face = true;
          if (!idx2.outOfBounds()) {
            auto cell2 = m_Cells.at(idx2);
            if (idx.level() == idx2.level() && cell.linear_index > cell2.linear_index) {
              use_face = false;
            }
          }
          if (use_face) {
            ++num_faces;
          }
        }
      }
    }
    return num_faces;
  }

  void extractFaces()
  {
    using namespace std;
    //
    m_Faces.clear();
    m_Faces.reserve(numPotentialFaces());
    amr_index_type neighbours[6];
    auto leaf_cells = getLeafCellIndices();
    for (auto idx : leaf_cells) {
      auto cell = m_Cells.at(idx);
      neighbours[0] = cellNeighbourIM(idx);
      neighbours[1] = cellNeighbourIP(idx);
      neighbours[2] = cellNeighbourJM(idx);
      neighbours[3] = cellNeighbourJP(idx);
      neighbours[4] = cellNeighbourKM(idx);
      neighbours[5] = cellNeighbourKP(idx);
      //
      for (int i = 0; i < 6; ++i) {
        auto idx2 = neighbours[i];
        if (idx2.valid() || idx2.outOfBounds()) {
          bool use_face = true;
          face_t face;
          face.nodes.resize(4);
          face.cell1 = cell.linear_index;          
          if (face.cell1 < 0) {
            throw std::runtime_error("AMRMesh: extractFaces: invalid cell mapping");
          }
          if (idx2.outOfBounds()) {
            face.cell2 = -1;
            if (cell.is_outside) {
              use_face = false;
            }
          } else {
            auto cell2 = m_Cells.at(idx2);
            if (cell.is_outside != cell2.is_outside) {
              face.is_boundary = true;
            }
            face.cell2 = cell2.linear_index;
            if (face.cell2 < 0) {
              throw std::runtime_error("AMRMesh: extractFaces: invalid cell mapping");
            }
            if (idx.level() == idx2.level() && cell.linear_index > cell2.linear_index) {
              use_face = false;
            }
            if (cell.is_outside && cell2.is_outside) {
              use_face = false;
            }
          }
          if (use_face) {
            switch (i) {
              case 0:
                face.nodes[0] = idx;
                face.nodes[1] = idx.incrementedK();
                face.nodes[2] = idx.incrementedJ().incrementedK();
                face.nodes[3] = idx.incrementedJ();
                break;
              case 1:
                face.nodes[0] = idx.incrementedI();
                face.nodes[1] = idx.incrementedI().incrementedJ();
                face.nodes[2] = idx.incrementedI().incrementedJ().incrementedK();
                face.nodes[3] = idx.incrementedI().incrementedK();
                break;
              case 2:
                face.nodes[0] = idx;
                face.nodes[1] = idx.incrementedI();
                face.nodes[2] = idx.incrementedI().incrementedK();
                face.nodes[3] = idx.incrementedK();
                break;
              case 3:
                face.nodes[0] = idx.incrementedJ();
                face.nodes[1] = idx.incrementedK().incrementedJ();
                face.nodes[2] = idx.incrementedI().incrementedJ().incrementedK();
                face.nodes[3] = idx.incrementedI().incrementedJ();
                break;
              case 4:
                face.nodes[0] = idx;
                face.nodes[1] = idx.incrementedJ();
                face.nodes[2] = idx.incrementedI().incrementedJ();
                face.nodes[3] = idx.incrementedI();
                break;
              case 5:
                face.nodes[0] = idx.incrementedK();
                face.nodes[1] = idx.incrementedI().incrementedK();
                face.nodes[2] = idx.incrementedI().incrementedJ().incrementedK();
                face.nodes[3] = idx.incrementedJ().incrementedK();
                break;
            }
            if (face.cell1 > face.cell2 && face.cell2 >= 0) {
              face.reverse();
            }
            completeFace(face);
            m_Faces.push_back(face);
          }
        }
      }
    }
    m_Faces.shrink_to_fit();
    //
    // vector<face_t> field_faces;
    // field_faces.reserve(m_Faces.size());
    // for (auto it : m_Faces) {
    //   if (it.cell2 >= 0) {
    //     field_faces.push_back(it);
    //   }
    // }
    // vector<face_t> boundary_faces;
    // boundary_faces.reserve(m_Faces.size() - field_faces.size());
    // for (auto it : m_Faces) {
    //   if (it.cell2 < 0) {
    //     boundary_faces.push_back(it);
    //   }
    // }
    // //
    // size_t i = 0;
    // for (const auto& face : field_faces) {
    //   m_Faces[i] = face;
    //   ++i;
    // }
    // for (const auto& face : boundary_faces) {
    //   m_Faces[i] = face;
    //   ++i;
    // }
    //
    // sort m_Faces
    //
    sort(m_Faces.begin(), m_Faces.end());
  }



public:

  AMRMesh(index_type sizeI, index_type sizeJ, index_type sizeK, const vector_type& x1, const vector_type& x2) : m_SizeI(sizeI), m_SizeJ(sizeJ), m_SizeK(sizeK), m_X1(x1), m_X2(x2)
  {
    for (size_t i = 0; i < m_SizeI; ++i) {
      for (size_t j = 0; j < m_SizeJ; ++j) {
        for (size_t k = 0; k < m_SizeK; ++k) {
          amr_index_type idx(i, j, k, 0);
          addCell(idx);
        }
      }
    }
  }

  void refineCell(const amr_index_type& idx)
  {
    auto it = m_Cells.find(idx);
    if (it == m_Cells.end()) {
      throw std::runtime_error("AMRMesh: cell not found");
    }
    amr_index_type idx0 = it->first.incrementedLevel();
    amr_index_type idx1 = idx0.incrementedI();
    amr_index_type idx2 = idx1.incrementedJ();
    amr_index_type idx3 = idx0.incrementedJ();
    amr_index_type idx4 = idx0.incrementedK();
    amr_index_type idx5 = idx1.incrementedK();
    amr_index_type idx6 = idx2.incrementedK();
    amr_index_type idx7 = idx3.incrementedK();
    addCell(idx0);
    addCell(idx1);
    addCell(idx2);
    addCell(idx3);
    addCell(idx4);
    addCell(idx5);
    addCell(idx6);
    addCell(idx7);
    m_Cells[idx].first_child = idx0;
  }

  void refineCell(index_type i, index_type j, index_type k, index_type level)
  {
    refineCell(amr_index_type(i, j, k, level));
  }

  index_type nodeLinearIndex(amr_index_type idx)  const
  {
    auto idx_it = m_Nodes.find(idx.reduced());
    if (idx_it == m_Nodes.end()) {
      throw std::runtime_error("AMRMesh: node not found");
    }
    return idx_it->second.linear_index;
  }

  size_t numCells() const
  {
    return m_Cells.size();
  }

  size_t numNodes() const
  {
    return m_Nodes.size();
  }

  void markCellAsOutside(const amr_index_type& idx)
  {
    auto it = m_Cells.find(idx);
    if (it == m_Cells.end()) {
      throw std::runtime_error("AMRMesh: cell not found");
    }
    it->second.is_outside = true;
  }

  bool cellMarkedAsOutside(const amr_index_type& idx) const
  {
    auto it = m_Cells.find(idx);
    if (it == m_Cells.end()) {
      throw std::runtime_error("AMRMesh: cell not found");
    }
    return it->second.is_outside;
  }

  std::vector<amr_index_type> nodesInBetweenI(const amr_index_type& idx0, const amr_index_type& idx1) const
  {
    if ((idx0.i() == idx1.i()) || (idx0.j() != idx1.j()) || (idx0.k() != idx1.k())) {
      throw std::runtime_error("AMRMesh: nodesInBetween: invalid index");
    }
    // std::vector<amr_index_type> nodes;
    // auto idx = idx0.indexOnLevel(m_MaxLevel).incrementedI();
    // while (idx.reduced() != idx1.reduced()) {
    //   if (m_Nodes.find(idx.reduced()) != m_Nodes.end()) { // if node idx exists
    //     nodes.push_back(idx);
    //   }
    //   idx.incrI();
    // }
    // return nodes;
    return nodesInBetween(idx0, idx1);
  }

  std::vector<amr_index_type> nodesInBetweenJ(const amr_index_type& idx0, const amr_index_type& idx1) const
  {
    if ((idx0.j() == idx1.j()) || (idx0.i() != idx1.i()) || (idx0.k() != idx1.k())) {
      throw std::runtime_error("AMRMesh: nodesInBetween: invalid index");
    }
    // std::vector<amr_index_type> nodes;
    // auto idx = idx0.indexOnLevel(m_MaxLevel).incrementedJ();
    // while (idx.reduced() != idx1.reduced()) {
    //   if (m_Nodes.find(idx.reduced()) != m_Nodes.end()) { // if node idx exists
    //     nodes.push_back(idx);
    //   }
    //   idx.incrJ();
    // }
    // return nodes;
    return nodesInBetween(idx0, idx1);
  }

  std::vector<amr_index_type> nodesInBetweenK(const amr_index_type& idx0, const amr_index_type& idx1) const
  {
    if ((idx0.k() == idx1.k()) || (idx0.i() != idx1.i()) || (idx0.j() != idx1.j())) {
      throw std::runtime_error("AMRMesh: nodesInBetween: invalid index");
    }
    // std::vector<amr_index_type> nodes;
    // auto idx = idx0.indexOnLevel(m_MaxLevel).incrementedK();
    // while (idx.reduced() != idx1.reduced()) {
    //   if (m_Nodes.find(idx.reduced()) != m_Nodes.end()) { // if node idx exists
    //     nodes.push_back(idx);
    //   }
    //   idx.incrK();
    // }
    // return nodes;
    return nodesInBetween(idx0, idx1);
  }

  std::vector<amr_index_type> nodesInBetween(amr_index_type idx0, amr_index_type idx1) const
  {
    using namespace std;
    //
    int di = idx1.i() - idx0.i();
    int dj = idx1.j() - idx0.j();
    int dk = idx1.k() - idx0.k();
    if (abs(di) + abs(dj) + abs(dk) != 1) {
      throw std::runtime_error("AMRMesh: nodesInBetween: nodes need to be aligned in one direction");
    }
    //
    vector<amr_index_type> nodes;
    idx0 = idx0.indexOnLevel(m_MaxLevel);
    amr_index_type idx(idx0.i() + di, idx0.j() + dj, idx0.k() + dk, idx0.level());
    while (idx.reduced() != idx1.reduced()) {
      if (m_Nodes.find(idx.reduced()) != m_Nodes.end()) { // if node idx exists
        nodes.push_back(idx);
      }
      idx = amr_index_type(idx.i() + di, idx.j() + dj, idx.k() + dk, idx.level());
    }
    return nodes;
  }

  amr_index_type cellNeighbourIP(const amr_index_type& idx) const
  {
    if (idx.i() == sizeIOnLevel(idx.level()) - 1) {
      return amr_index_type::out_of_bounds();
    }
    amr_index_type I = cellNeighbour(idx, idx.incrementedI());
    return I;
  }

  amr_index_type cellNeighbourJP(const amr_index_type& idx) const
  {
    if (idx.j() == sizeJOnLevel(idx.level()) - 1) {
      return amr_index_type::out_of_bounds();
    }
    return cellNeighbour(idx, idx.incrementedJ());
  }

  amr_index_type cellNeighbourKP(const amr_index_type& idx) const
  {
    if (idx.k() == sizeKOnLevel(idx.level()) - 1) {
      return amr_index_type::out_of_bounds();
    }
    return cellNeighbour(idx, idx.incrementedK());
  }

  amr_index_type cellNeighbourIM(const amr_index_type& idx) const
  {
    if (idx.i() == 0) {
      return amr_index_type::out_of_bounds();
    }
    return cellNeighbour(idx, idx.decrementedI());
  }

  amr_index_type cellNeighbourJM(const amr_index_type& idx) const
  {
    if (idx.j() == 0) {
      return amr_index_type::out_of_bounds();
    }
    return cellNeighbour(idx, idx.decrementedJ());
  }

  amr_index_type cellNeighbourKM(const amr_index_type& idx) const
  {
    if (idx.k() == 0) {
      return amr_index_type::out_of_bounds();
    }
    return cellNeighbour(idx, idx.decrementedK());
  }

  std::vector<int> cellMapping() const
  {
    std::vector<amr_index_type> i2idx(m_Cells.size());
    for (auto it : m_Cells) {
      auto idx = it.first;
      auto cell = it.second;
      i2idx[cell.linear_index] = idx;
    }
    std::vector<int> mapping(m_Cells.size(), -1);
    size_t j = 0;
    for (size_t i = 0; i < i2idx.size(); ++i) {
      if (m_Cells.find(i2idx[i]) == m_Cells.end()) {
        throw std::runtime_error("AMRMesh: cellMapping: invalid index");
      }
      auto cell = m_Cells.at(i2idx[i]);
      if (!cell.first_child.valid()) {
        mapping[i] = j;
        ++j;
      }
    }
    return mapping;
  }

  vector_type getNodeCoordinates(const amr_index_type& idx) const
  {
    auto wi = typename vector_type::value_type(idx.i());
    auto wj = typename vector_type::value_type(idx.j());
    auto wk = typename vector_type::value_type(idx.k());
    vector_type point;
    point[0] = m_X1[0] + (m_X2[0] - m_X1[0]) * (wi / sizeIOnLevel(idx.level()));
    point[1] = m_X1[1] + (m_X2[1] - m_X1[1]) * (wj / sizeJOnLevel(idx.level()));
    point[2] = m_X1[2] + (m_X2[2] - m_X1[2]) * (wk / sizeKOnLevel(idx.level()));
    return point;
  }

  std::vector<vector_type> extractPoints() const
  {
    using namespace std;
    //
    vector<vector_type> points(m_Nodes.size());
    for (auto it : m_Nodes) {
      auto idx  = it.first;
      auto node = it.second;
      auto i    = node.linear_index;
      points[i] = getNodeCoordinates(idx);
    }
    return points;
  }

  std::vector<amr_index_type> getNodeIndices() const
  {
    std::vector<amr_index_type> indices;
    indices.reserve(m_Nodes.size());
    for (auto it : m_Nodes) {
      indices.push_back(it.first);
    }
    return indices;
  }

  std::vector<amr_index_type> getLeafCellIndices()
  {
    auto cells_backup = m_Cells;
    //
    // set all linear indices to -1
    //
    for (auto it : m_Cells) {
      auto idx = it.first;
      auto cell = it.second;
      cell.linear_index = -1;
      m_Cells[idx] = cell;
    }
    //
    std::vector<amr_index_type> indices;
    indices.reserve(m_Cells.size());
    //
    // add all leaf cells which are not outside 
    // and set the linear index to the index in the vector 
    //
    for (auto it : m_Cells) {
      auto idx  = it.first;
      auto cell = it.second;
      if (!cell.first_child.valid() && !cell.is_outside) {
        indices.push_back(idx);
        cell.linear_index = indices.size() - 1;
        m_Cells[idx] = cell;
      }
    }
    //
    // add all leaf cells which are outside
    // and set the linear index to the index in the vector
    //
    for (auto it : m_Cells) {
      auto idx  = it.first;
      auto cell = it.second;
      if (!cell.first_child.valid() && cell.is_outside) {
        indices.push_back(idx);
        cell.linear_index = indices.size() - 1;
        m_Cells[idx] = cell;
      }
    }
    //
    // set all the other linear indices to unique values
    //
    index_type linear_index = indices.size();
    for (auto it : m_Cells) {
      auto idx = it.first;
      auto cell = it.second;
      if (cell.linear_index < 0) {
        cell.linear_index = linear_index;
        m_Cells[idx] = cell;
        ++linear_index;
      }
    }
    //
    indices.shrink_to_fit();
    return indices;
  }

  vector_type cellCoordinates(const amr_index_type& idx) const
  {
    auto wi = typename vector_type::value_type(idx.i()) + 0.5;
    auto wj = typename vector_type::value_type(idx.j()) + 0.5;
    auto wk = typename vector_type::value_type(idx.k()) + 0.5;
    vector_type point;
    point[0] = m_X1[0] + (m_X2[0] - m_X1[0]) * (wi / sizeIOnLevel(idx.level()));
    point[1] = m_X1[1] + (m_X2[1] - m_X1[1]) * (wj / sizeJOnLevel(idx.level()));
    point[2] = m_X1[2] + (m_X2[2] - m_X1[2]) * (wk / sizeKOnLevel(idx.level()));
    return point;
  }

  std::vector<vector_type> getNodeCoordinatesOfCell(amr_index_type cell_idx) const
  {
    std::vector<vector_type> coords;
    coords.reserve(8);
    auto node0 = cell_idx;
    auto node1 = node0.incrementedI();
    auto node2 = node1.incrementedJ();
    auto node3 = node0.incrementedJ();
    auto node4 = node0.incrementedK();
    auto node5 = node1.incrementedK();
    auto node6 = node2.incrementedK();
    auto node7 = node3.incrementedK();
    coords.push_back(getNodeCoordinates(node0));
    coords.push_back(getNodeCoordinates(node1));
    coords.push_back(getNodeCoordinates(node2));
    coords.push_back(getNodeCoordinates(node3));
    coords.push_back(getNodeCoordinates(node4));
    coords.push_back(getNodeCoordinates(node5));
    coords.push_back(getNodeCoordinates(node6));
    coords.push_back(getNodeCoordinates(node7));
    return coords;
  }

  void writeFoamPointsFile(const std::string& file_name, std::vector<vector_type> X=std::vector<vector_type>()) const
  {
    using namespace std;
    //
    if (X.size() == 0) {
      X = extractPoints();
    }
    cout << "writing " << X.size() << " points to file " << file_name << endl;
    //
    ofstream file(file_name);
    file << "FoamFile\n";
    file << "{\n";
    file << "    version     2.0;\n";
    file << "    format      ascii;\n";
    file << "    class       vectorField;\n";
    file << "    location    \"constant/polyMesh\";\n";
    file << "    object      points;\n";
    file << "}\n\n";
    file << X.size() << "\n";
    file << "(\n";
    for (const auto &point : X) {
      file << "    (" << point[0] << " " << point[1] << " " << point[2] << ")\n";
    }
    file << ")\n";
  }

  void writeFoamFacesFile(const std::string& file_name) const
  {
    using namespace std;
    //
    cout << "writing " << m_Faces.size() << " faces to file " << file_name << endl;
    //
    ofstream file(file_name);
    file << "FoamFile\n";
    file << "{\n";
    file << "    version     2.0;\n";
    file << "    format      ascii;\n";
    file << "    class       faceList;\n";
    file << "    location    \"constant/polyMesh\";\n";
    file << "    object      faces;\n";
    file << "}\n\n";
    file << m_Faces.size() << "\n";
    file << "(\n";
    for (const auto &face : m_Faces) {
      file << face.nodes.size() << "(";
      for (const auto &node : face.nodes) {
        auto n = m_Nodes.at(node.reduced()).linear_index;
        file << " " << n;
      }
      file << ")\n";
    }
    file << ")\n";
  }

  void writeFoamOwnerFile(const std::string& file_name) const
  {
    using namespace std;
    //
    ofstream file(file_name);
    file << "FoamFile\n";
    file << "{\n";
    file << "    version     2.0;\n";
    file << "    format      ascii;\n";
    file << "    class       labelList;\n";
    file << "    location    \"constant/polyMesh\";\n";
    file << "    object      owner;\n";
    file << "}\n\n";
    file << m_Faces.size() << "\n";
    file << "(\n";
    for (const auto &face : m_Faces) {
      file << face.cell1 << "\n";
    }
    file << ")\n";
  }

  void writeFoamNeighbourFile(const std::string& file_name) const
  {
    using namespace std;
    //
    ofstream file(file_name);
    file << "FoamFile\n";
    file << "{\n";
    file << "    version     2.0;\n";
    file << "    format      ascii;\n";
    file << "    class       labelList;\n";
    file << "    location    \"constant/polyMesh\";\n";
    file << "    object      neighbour;\n";
    file << "}\n\n";
    int N = 0;
    for (const auto &face : m_Faces) {
      if (face.cell2 >= 0) {
        ++N;
      }
    }
    file << N << "\n";  
    file << "(\n";
    for (const auto &face : m_Faces) {
      if (face.cell2 >= 0) {
        file << face.cell2 << "\n";
      }
    }
    file << ")\n";
  }

  void writeFoamBoundaryFile(const std::string& file_name) const
  {
    using namespace std;
    //
    int n_faces    = 0;
    int start_face = 0;
    for (const auto &face : m_Faces) {
      if (face.cell2 >= 0) {
        ++start_face;
      } else {
        ++n_faces;
      }
    }
    ofstream file(file_name);
    file << "FoamFile\n";
    file << "{\n";
    file << "    version     2.0;\n";
    file << "    format      ascii;\n";
    file << "    class       polyBoundaryMesh;\n";
    file << "    location    \"constant/polyMesh\";\n";
    file << "    object      boundary;\n";
    file << "}\n\n";
    file << "1\n";
    file << "(\n";
    file << "    all\n";
    file << "    {\n";
    file << "        type patch;\n";
    file << "        nFaces    " << n_faces << ";\n";
    file << "        startFace " << start_face << ";\n";
    file << "    }\n";
    file << ")\n";
  }
  
  void writeFoamMesh(const std::string& path, std::vector<vector_type> X=std::vector<vector_type>())
  {
    using namespace std;
    //
    extractFaces();
    //
    // sort faces to match OpenFOAM requirements
    //
    std::vector<face_t> field_faces;
    field_faces.reserve(m_Faces.size());
    for (auto it : m_Faces) {
      if (it.cell2 >= 0) {
        field_faces.push_back(it);
      }
    }
    {
      std::vector<face_t> boundary_faces;
      boundary_faces.reserve(m_Faces.size() - field_faces.size());
      for (auto it : m_Faces) {
        if (it.cell2 < 0) {
          boundary_faces.push_back(it);
        }
      }
      //
      size_t i = 0;
      for (const auto& face : field_faces) {
        m_Faces[i] = face;
        ++i;
      }
      for (const auto& face : boundary_faces) {
        m_Faces[i] = face;
        ++i;
      }
    }
    //
    writeFoamPointsFile(path + "/points", X);
    writeFoamFacesFile(path + "/faces");
    writeFoamOwnerFile(path + "/owner");
    writeFoamNeighbourFile(path + "/neighbour");
    writeFoamBoundaryFile(path + "/boundary");
  }

  std::vector<face_t>& faces()
  {
    if (m_Faces.size() == 0) {
      extractFaces();
    }
    return m_Faces;
  }

  std::vector<amr_index_type> cellNeighbourIndices(amr_index_type idx) const
  {
    std::vector<amr_index_type> indices, neigh(6);
    neigh[0] = cellNeighbourIM(idx);
    neigh[1] = cellNeighbourIP(idx);
    neigh[2] = cellNeighbourJM(idx);
    neigh[3] = cellNeighbourJP(idx);
    neigh[4] = cellNeighbourKM(idx);
    neigh[5] = cellNeighbourKP(idx);
    indices.reserve(6);
    //
    for (auto _idx : neigh) {
      if (_idx.valid() && !_idx.outOfBounds()) {
        indices.push_back(_idx);
      }
    }
    //
    return indices;
  }

  void ensureSmoothTransition(int num_layers=1)
  {
    using namespace std;
    //    
    bool done = false;
    //
    // count leaf cells
    //
    size_t num_leaf_cells = 0;
    for (auto it : m_Cells) {
      auto idx  = it.first;
      auto cell = it.second;
      if (!cell.first_child.valid()) {
        ++num_leaf_cells;
      }
    }
    cout << "AMRMesh: ensureSmoothTransition: starting with " << num_leaf_cells << " cells" << endl;
    amr_index_type neigh[6];
    while (!done) {
      done = true;
      int N = 0;
      //
      // loop over all cells and check if they have a neighbour with a level that is more than one level higher
      //
      for (auto it : m_Cells) {
        auto idx  = it.first;
        auto cell = it.second;
        if (cell.first_child.valid()) {
          continue;
        }
        neigh[0] = cellNeighbourIM(idx);
        neigh[1] = cellNeighbourIP(idx);
        neigh[2] = cellNeighbourJM(idx);
        neigh[3] = cellNeighbourJP(idx);
        neigh[4] = cellNeighbourKM(idx);
        neigh[5] = cellNeighbourKP(idx);
        for (int i = 0; i < 6; ++i) {
          auto idx2 = neigh[i];
          if (idx2.valid() && idx2.level() < idx.level() - 1) {
            refineCell(idx2);
            done = false;
            ++N;
          }
        }
      }
      cout << "AMRMesh: ensureSmoothTransition: added " << N << " cells" << endl;
    }
    //
    cout << "AMRMesh: ensureSmoothTransition: done" << endl;
  }

  void ensureSmoothTransition2()
  {
    using namespace std;
    //    
    bool done = false;
    //
    // count leaf cells
    //
    size_t num_leaf_cells = 0;
    for (auto it : m_Cells) {
      auto idx  = it.first;
      auto cell = it.second;
      if (!cell.first_child.valid()) {
        ++num_leaf_cells;
      }
    }
    cout << "AMRMesh: ensureSmoothTransition: starting with " << num_leaf_cells << " cells" << endl;
    //
    // mark all nodes which belong to a cell of the maximal refinement level
    //
    for (int level = m_MaxLevel; level >= 0; --level) {
      unordered_map<amr_index_type,bool> nodes_on_level;
      while (!done) {
        done = true;
        int N = 0;
      }
    }
  }

  bool nodeOnXBoundary(const amr_index_type& idx) const
  {
    return idx.i() == 0 || idx.i() == sizeIOnLevel(idx.level());
  }

  bool nodeOnYBoundary(const amr_index_type& idx) const
  {
    return idx.j() == 0 || idx.j() == sizeJOnLevel(idx.level());
  }

  bool nodeOnZBoundary(const amr_index_type& idx) const
  {
    return idx.k() == 0 || idx.k() == sizeKOnLevel(idx.level());
  }

  bool nodeOnMinusXBoundary(const amr_index_type& idx) const
  {
    return idx.i() == 0;
  }

  bool nodeOnPlusXBoundary(const amr_index_type& idx) const
  {
    return idx.i() == sizeIOnLevel(idx.level());
  }

  bool nodeOnMinusYBoundary(const amr_index_type& idx) const
  {
    return idx.j() == 0;
  }

  bool nodeOnPlusYBoundary(const amr_index_type& idx) const
  {
    return idx.j() == sizeJOnLevel(idx.level());
  }

  bool nodeOnMinusZBoundary(const amr_index_type& idx) const
  {
    return idx.k() == 0;
  }

  bool nodeOnPlusZBoundary(const amr_index_type& idx) const
  {
    return idx.k() == sizeKOnLevel(idx.level());
  }

};

} // namespace edl


namespace std {
  template <typename TIndex>
  struct hash<edl::AMRIndex<TIndex>>
  {
    size_t operator()(const edl::AMRIndex<TIndex>& node) const
    {
      // static const uint8_t p1 = 31;
      // static const uint8_t p2 = 37;
      // static const uint8_t p3 = 41;
      // static const uint8_t p4 = 43;
      static const uint64_t p1 = 1;
      static const uint64_t p2 = 65536;
      static const uint64_t p3 = 4294967296;
      static const uint64_t p4 = 281474976710656;
      return node.i()*p1 + node.j()*p2 + node.k()*p3 + node.level()*p4;
    }
  };
}


// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

TEST_CASE("AMRIndex") 
{
  using namespace std;
  using namespace edl;
  //
  typedef AMRIndex<int32_t> idx_t;
  idx_t n1(8, 16, 24, 3);
  idx_t n2(1, 2, 3, 0);
  CHECK(n1.reduced() == n2);
  CHECK(idx_t(8, 16, 25,3).reduced().level() == 3);
  auto n3 = n1 + n2;
  auto n4 = n1 - n2;
  CHECK(n3 == 2*n1);
  CHECK(n3/2 == n1);
  CHECK(n4.reduced() == idx_t(0,0,0));
  CHECK(n1/4 == idx_t(2,4,6,3));
  CHECK(n2/4 == idx_t(1,2,3,2));
  //
  idx_t q1(10,10,3,0);
  idx_t q2(11,10,3,0);
  idx_t q3(11,11,3,0);
  idx_t q4(10,11,3,0);
  idx_t qc = (q1 + q2 + q3 + q4) / 4;
  CHECK(qc == idx_t(21,21,6,1));
}

TEST_CASE("AMRMesh_basics") 
{
  using namespace std;
  using namespace edl;
  //
  typedef MathVector<StaticVector<float,3> > vec3_t;
  typedef AMRMesh<int32_t, uint16_t, vec3_t> mesh_t;
  //
  int N = 10;
  mesh_t mesh(N, N, N, vec3_t(0,0,0), vec3_t(1,1,1));
  CHECK(mesh.numCells() == N*N*N);
  CHECK(mesh.numNodes() == (N+1)*(N+1)*(N+1));  
  //
  mesh_t::amr_index_type idx1(2,2,2,0);
  mesh.refineCell(idx1);
  CHECK(mesh.numCells() == N*N*N + 8);
  CHECK(mesh.numNodes() == (N+1)*(N+1)*(N+1) + 12 + 6 + 1);
  auto idx2I = idx1.incrementedI();
  auto nodes_in_between1 = mesh.nodesInBetweenI(idx1, idx2I);
  CHECK(nodes_in_between1.size() == 1);
  auto idx3 = idx1.indexOnLevel(1);
  mesh.refineCell(idx3);
  auto nodes_in_between2 = mesh.nodesInBetweenI(idx1, idx2I);
  CHECK(nodes_in_between2.size() == 2);
  auto idx2J = idx1.incrementedJ();
  CHECK(mesh.nodesInBetweenJ(idx1, idx2J).size() == 2);
  auto idx2K = idx1.incrementedK();
  CHECK(mesh.nodesInBetweenK(idx1, idx2K).size() == 2);
}

TEST_CASE("AMRMesh_cell_neighbours") 
{
  using namespace std;
  using namespace edl;
  //
  typedef MathVector<StaticVector<float,3> > vec3_t;
  typedef AMRMesh<int32_t, uint16_t, vec3_t> mesh_t;
  typedef mesh_t::amr_index_type idx_t;
  //
  int N = 10;
  mesh_t mesh(N, N, N, vec3_t(0,0,0), vec3_t(1,1,1));
  //
  idx_t I1(1,1,1);
  idx_t I2 = I1.indexOnLevel(1);
  idx_t I3 = I2.indexOnLevel(2);
  mesh.refineCell(I1);
  mesh.refineCell(I2);
  //
  auto I3_nip = mesh.cellNeighbourIP(I3);
  CHECK(I3_nip == idx_t(5,4,4,2));
  auto I3_njp = mesh.cellNeighbourJP(I3);
  CHECK(I3_njp == idx_t(4,5,4,2));
  auto I3_nkp = mesh.cellNeighbourKP(I3);
  CHECK(I3_nkp == idx_t(4,4,5,2));
  auto tmp = mesh.cellNeighbourIM(I3_nip);
  CHECK(mesh.cellNeighbourIM(I3_nip) == I3);
  CHECK(mesh.cellNeighbourJM(I3_njp) == I3);
  CHECK(mesh.cellNeighbourKM(I3_nkp) == I3);
  auto dbg = mesh.cellNeighbourIM(I3);
  CHECK(mesh.cellNeighbourIM(I3) == idx_t(0,1,1,0));
  //
  auto I4 = I3_nip;
  {
    auto tmp = mesh.cellNeighbourIP(I4);
    CHECK(tmp == idx_t(3,2,2,1));
    CHECK(mesh.cellNeighbourIM(tmp) == idx_t::invalid());
  }
  auto I5 = I3_njp;
  {
    auto tmp = mesh.cellNeighbourJP(I5);
    CHECK(tmp == idx_t(2,3,2,1));
    CHECK(mesh.cellNeighbourJM(tmp) == idx_t::invalid());
  }
  auto I6 = I3_nkp;
  {
    auto tmp = mesh.cellNeighbourKP(I6);
    CHECK(tmp == idx_t(2,2,3,1));
    CHECK(mesh.cellNeighbourKM(tmp) == idx_t::invalid());
  }
  //
  idx_t I7(0,0,0);
  {
    auto tmp = mesh.cellNeighbourIM(I7);
    CHECK(tmp == idx_t::out_of_bounds());
  }
  {
    auto tmp = mesh.cellNeighbourJM(I7);
    CHECK(tmp == idx_t::out_of_bounds());
  }
  {
    auto tmp = mesh.cellNeighbourKM(I7);
    CHECK(tmp == idx_t::out_of_bounds());
  }
  idx_t I8(N-1,N-1,N-1);
  {
    auto tmp = mesh.cellNeighbourIP(I8);
    CHECK(tmp == idx_t::out_of_bounds());
  }
  {
    auto tmp = mesh.cellNeighbourJP(I8);
    CHECK(tmp == idx_t::out_of_bounds());
  }
  {
    auto tmp = mesh.cellNeighbourKP(I8);
    CHECK(tmp == idx_t::out_of_bounds());
  }
}

TEST_CASE("AMRMesh_extract_faces")
{
  using namespace std;
  using namespace edl;
  //
  typedef MathVector<StaticVector<float,3> > vec3_t;
  typedef AMRMesh<int32_t, uint16_t, vec3_t> mesh_t;
  typedef mesh_t::amr_index_type idx_t;
  //
  int N = 2;
  mesh_t mesh(N, N, N, vec3_t(0,0,0), vec3_t(1,1,1));
  mesh.refineCell(0,0,0,0);
  auto faces = mesh.faces();
  for (int i = 0; i < faces.size(); ++i) {
    cout << "face " << i << endl;
    cout << "  - cell1 : " << faces[i].cell1 << endl;
    cout << "  - cell2 : " << faces[i].cell2 << endl;
    cout << "  - nodes  : ";
    for (int j = 0; j < 4; ++j) {
      cout << faces[i].nodes[j] << " ";
    }
    cout << endl;
  }
}

TEST_CASE("AMRMesh_sphere_example")
{
  using namespace std;
  using namespace edl;
  //
  typedef MathVector<StaticVector<float,3> > vec3_t;
  typedef AMRMesh<int32_t, uint16_t, vec3_t> mesh_t;
  typedef mesh_t::amr_index_type idx_t;
  //
  int N = 10;
  mesh_t mesh(N, N, N, vec3_t(-2,-2,-2), vec3_t(2,2,2));
  //
  // refine cells which cross the surface of a sphere at (0,0,0) with radius 1
  //
  for (int iter = 0; iter < 3; ++iter) {
    auto cells = mesh.getLeafCellIndices();
    for (auto idx : cells) {
      auto coords = mesh.getNodeCoordinatesOfCell(idx);
      bool pos    = false;
      bool neg    = false;
      for (auto& coord : coords) {
        if (coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2] < 1) {
          neg = true;
        } else {
          pos = true;        
        }
      }
      if (pos && neg) {
        mesh.refineCell(idx);
      }
    }
  }
  mesh.ensureSmoothTransition();
}

TEST_CASE("AMRMesh_cell_neighbours_for_layers") 
{
  using namespace std;
  using namespace edl;
  //
  typedef MathVector<StaticVector<float,3> > vec3_t;
  typedef AMRMesh<int32_t, uint16_t, vec3_t> mesh_t;
  typedef mesh_t::amr_index_type idx_t;
  //
  int N = 5;
  mesh_t mesh(N, N, N, vec3_t(0,0,0), vec3_t(1,1,1));
  //
  mesh.refineCell(idx_t(2,2,2,0));
  //
  system("mkdir -p constant/polyMesh");
  mesh.writeFoamMesh("constant/polyMesh");
  ofstream file("test.foam");
  file.close();
  //
  idx_t I_212_0(2,1,2, 0);
  idx_t I_222_0(2,2,2, 0);
  idx_t I_444_1(4,4,4, 1);
  //
  CHECK(I_444_1.indexOnLevel(0) == I_222_0);
  //
  auto I_212_0_nyp = mesh.cellNeighbourJP(I_212_0);
  auto I_212_0_nym = mesh.cellNeighbourJM(I_212_0);
  auto I_222_0_nyp = mesh.cellNeighbourJP(I_222_0);
  auto I_222_0_nym = mesh.cellNeighbourJM(I_222_0);
  auto I_444_1_nyp = mesh.cellNeighbourJP(I_444_1);
  auto I_444_1_nym = mesh.cellNeighbourJM(I_444_1);
  //
  CHECK(!I_212_0_nyp.valid());
  CHECK(I_212_0_nym == idx_t(2,0,2,0));
  CHECK(I_222_0_nyp == idx_t(2,3,2,0));
  CHECK(I_222_0_nym == idx_t(2,1,2,0));
  CHECK(I_444_1_nyp == idx_t(4,5,4,1));
  CHECK(I_444_1_nym == idx_t(2,1,2,0));
  //
  mesh.refineCell(I_444_1);
  auto I_888_2 = I_444_1.indexOnLevel(2);
  auto I_888_2_nyp = mesh.cellNeighbourJP(I_888_2);
  auto I_888_2_nym = mesh.cellNeighbourJM(I_888_2);
  CHECK(I_888_2 == idx_t(8,8,8,2));
  CHECK(I_888_2_nyp == idx_t(8,9,8,2));
  CHECK(I_888_2_nym == idx_t(2,1,2,0));
}


#endif // AMRINDEX_H