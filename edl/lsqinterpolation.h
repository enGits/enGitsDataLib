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

#ifndef LSQINTERPOLATION_H
#define LSQINTERPOLATION_H

#include "edl/edl.h"
#include "edl/mathvector.h"
#include "edl/smallsquarematrix.h"

#include <vector>

namespace EDL_NAMESPACE
{
template <typename T, uint8_t DIM> class LsqInterpolation;
}

namespace EDL_NAMESPACE
{

template <typename VALUE_TYPE, uint8_t DIM>
class LsqInterpolation
{

protected: // data types

  typedef VALUE_TYPE                              value_t;
  typedef MathVector<StaticVector<value_t,DIM>>   geovec_t;
  typedef MathVector<StaticVector<value_t,DIM+1>> vector_t;
  typedef SmallSquareMatrix<value_t,DIM+1>        matrix_t;

  struct node_t
  {
    geovec_t x;
    value_t  weight;
    value_t  value;
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
      A[0][0] += node.weight;
      for (int ij = 1; ij <= DIM; ++ij) {
        A[ij][0] += node.weight*node.x[ij-1];
        A[0][ij] += node.weight*node.x[ij-1];
      }
      for (int i = 1; i <= DIM; ++i) {
        for (int j = 1; j <= DIM; ++j) {
          A[i][j] += node.weight*node.x[i-1]*node.x[j-1];
        }
      }
    }
    m_Matrix = A;
    m_InvMatrix = A.inverse();
    //
    m_MatrixPrepared = true;
  }

  void prepareVector()
  {
    if (!m_MatrixPrepared) {
      prepareMatrix();
    }
    for (int i = 0; i <= DIM; ++i) {
      m_B[i] = 0;
    }
    for (auto node : m_Nodes) {
      m_B[0] += node.weight*node.value;
      for (int i = 1; i <= DIM; ++i) {
        m_B[i] += node.weight*node.value*node.x[i-1];
      }
    }
    m_A = m_InvMatrix*m_B;
    m_VectorPrepared = true;
  }


public: // methods

  LsqInterpolation()
  {
    // pass
  }

  template<typename T>
  void addNode(const T& x, value_t w=1)
  {
    node_t node;
    node.weight = w;
    node.value  = 0;
    for (int i = 0; i < DIM; ++i) {
      node.x[i] = x[i];
    }
    m_Nodes.push_back(node);
    m_MatrixPrepared = false;
    m_VectorPrepared = false;
  }

  void setValue(int i, value_t value)
  {
    m_VectorPrepared = false;
    m_Nodes[i].value = value;
  }

  template <typename T>
  value_t interpolate(const T& x)
  {
    if (!m_VectorPrepared) {
      prepareVector();
    }
    value_t value = m_A[0];
    for (int i = 1; i <= DIM; ++i) {
      value += x[i-1]*m_A[i];
    }
    return value;
  }

  void reset()
  {
    m_Nodes.clear();
  }

  int size()
  {
    return m_Nodes.size();
  }

  vector_t coeffs()
  {
    if (!m_VectorPrepared) {
      prepareVector();
    }
    return m_A;
  }

  void write(std::ostream& s)
  {
    if (!m_VectorPrepared) {
      prepareVector();
    }
    for (int i = 0; i <= DIM; ++i) {
      for (int j = 0; j <= DIM; ++j) {
        s << m_Matrix[i][j] << ", ";
      }
      s << m_B[i] << "\n";
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

  template<typename T>
  std::vector<value_t> getMetrics(const T& x)
  {
    if (!m_MatrixPrepared) {
      prepareMatrix();
    }
    //
    std::vector<value_t> C(m_Nodes.size());
    //a_00 + a_01*x[i] + a_02*y[i] + a_03*z[i] + a_10*x + a_11*x*x[i] + a_12*x*y[i] + a_13*x*z[i] + a_20*y + a_21*x[i]*y + a_22*y*y[i] + a_23*y*z[i] + a_30*z + a_31*x[i]*z + a_32*y[i]*z + a_33*z*z[i]

    for (int i = 0; i < m_Nodes.size(); ++i) {
      node_t node = m_Nodes[i];
      C[i]  = m_InvMatrix[0][0];
      C[i] += m_InvMatrix[0][1]*node.x[0];
      C[i] += m_InvMatrix[0][2]*node.x[1];
      C[i] += m_InvMatrix[0][3]*node.x[2];
      C[i] += m_InvMatrix[1][0]*x[0];
      C[i] += m_InvMatrix[1][1]*x[0]*node.x[0];
      C[i] += m_InvMatrix[1][2]*x[0]*node.x[1];
      C[i] += m_InvMatrix[1][3]*x[0]*node.x[2];
      C[i] += m_InvMatrix[2][0]*x[1];
      C[i] += m_InvMatrix[2][1]*x[1]*node.x[0];
      C[i] += m_InvMatrix[2][2]*x[1]*node.x[1];
      C[i] += m_InvMatrix[2][3]*x[1]*node.x[2];
      C[i] += m_InvMatrix[3][0]*x[2];
      C[i] += m_InvMatrix[3][1]*x[2]*node.x[0];
      C[i] += m_InvMatrix[3][2]*x[2]*node.x[1];
      C[i] += m_InvMatrix[3][3]*x[2]*node.x[2];
      C[i] *= node.weight;
    }
    return C;
  }

};


} // EDL_NAMESPACE

#endif // TGRAPH_H


