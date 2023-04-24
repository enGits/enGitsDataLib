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

#ifndef LSQGRAD_H
#define LSQGRAD_H

#include "edl/edl.h"
#include "edl/mathvector.h"
#include "edl/smallsquarematrix.h"

#include <vector>

namespace EDL_NAMESPACE
{
template <typename T, uint8_t DIM> class LsqGrad;
}

namespace EDL_NAMESPACE
{

template <typename VALUE_TYPE, uint8_t DIM>
class LsqGrad
{

protected: // data types

  typedef VALUE_TYPE                            value_t;
  typedef MathVector<StaticVector<value_t,DIM>> vector_t;
  typedef SmallSquareMatrix<value_t,DIM>        matrix_t;

  struct node_t
  {
    vector_t x;
    value_t  weight;
  };


protected: // attributes

  std::vector<node_t> m_Nodes;
  matrix_t            m_Matrix;
  matrix_t            m_InvMatrix;
  vector_t            m_B;
  vector_t            m_A;
  bool                m_MatrixPrepared = false;
  bool                m_VectorPrepared = false;


protected: // methods

  void prepareMatrix()
  {
    matrix_t A;
    A.initAll(0);
    //
    for (auto node : m_Nodes) {
      for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) {
          A[i][j]+= node.weight*node.x[i]*node.x[j];
        }
      }
    }
    m_Matrix = A;
    m_InvMatrix = A.inverse();
    //
    m_MatrixPrepared = true;
  }


public: // methods

  template<typename T>
  void addNode(const T& x, value_t w=1)
  {
    node_t node;
    node.weight = w;
    for (int i = 0; i < DIM; ++i) {
      node.x[i] = x[i];
    }
    m_Nodes.push_back(node);
    m_MatrixPrepared = false;
    m_VectorPrepared = false;
  }

  void reset()
  {
    m_Nodes.clear();
  }

  int size()
  {
    return m_Nodes.size();
  }

  void write(std::ostream& s)
  {
    for (int i = 0; i <= DIM; ++i) {
      for (int j = 0; j <= DIM; ++j) {
        s << m_Matrix[i][j] << ", ";
      }
      s << "\n";
    }
  }

  void getFlatMatrix(value_t* flat_matrix)
  {
    if (!m_MatrixPrepared) {
      prepareMatrix();
    }
    int idx = 0;
    for (int i = 0; i <= DIM; ++i) {
      for (int j = 0; j <= DIM; ++j) {
        flat_matrix[idx] = m_InvMatrix[i][j];
        ++idx;
      }
    }
  }

  matrix_t getMatrix()
  {
    return m_Matrix;
  }

  matrix_t getInvMatrix()
  {
    return m_InvMatrix;
  }

  std::vector<vector_t> getMetrics()
  {
    if (!m_MatrixPrepared) {
      prepareMatrix();
    }
  
    std::vector<std::vector<value_t>> C( m_Nodes.size(), std::vector<value_t>(3));
    for (int i = 0; i < m_Nodes.size(); ++i) {
      node_t node = m_Nodes[i];
      for (int j = 0; j < DIM; ++j) {
        for (int k = 0; k < DIM; ++k) {
          C[i][k] += (m_InvMatrix[k][j] * node.x[j]) * node.weight;
        }
      }
    }
    return C;
  }

};


} // EDL_NAMESPACE

#endif // LSQGRAD_H


