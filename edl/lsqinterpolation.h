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
template <typename T> class LsqInterpolation3D;
}

namespace EDL_NAMESPACE
{

template <typename VALUE_TYPE>
class LsqInterpolation3D
{

protected: // data types

  typedef VALUE_TYPE                          value_t;
  typedef MathVector<StaticVector<value_t,3>> geovec_t;
  typedef MathVector<StaticVector<value_t,4>> vector_t;
  typedef SmallSquareMatrix<value_t,4>        matrix_t;


protected: // attributes

  std::vector<value_t> m_Alpha;


public: // methods


  template<typename C>
  void setNodes(geovec_t x0, C& X, value_t distance_exponent=0.0)
  {
    value_t min_dist = std::numeric_limits<value_t>::max();
    value_t max_dist = 0;
    for (auto x : X) {
      min_dist = std::min(min_dist, (x - x0).abs());
      max_dist = std::max(max_dist, (x - x0).abs());
    }
    min_dist = std::max(min_dist, max_dist/100);
    m_Alpha.resize(X.size());
    matrix_t A;
    std::vector<value_t> weights;
    weights.reserve(X.size());
    A.initAll(0);
    for (auto x : X) {
      value_t  w  = std::pow((x - x0).abs()/min_dist, distance_exponent);
      vector_t DX = x - x0;
      value_t  dx = DX[0];
      value_t  dy = DX[1];
      value_t  dz = DX[2];
      //
      weights.push_back(w);
      //
      A[0][0] += 2*w;
      A[0][1] += 2*dx*w;
      A[0][2] += 2*dy*w;
      A[0][3] += 2*dz*w;
      A[1][0] += 2*dx*w;
      A[1][1] += 2*dx*dx*w;
      A[1][2] += 2*dx*dy*w;
      A[1][3] += 2*dx*dz*w;
      A[2][0] += 2*dy*w;
      A[2][1] += 2*dx*dy*w;
      A[2][2] += 2*dy*dy*w;
      A[2][3] += 2*dy*dz*w;
      A[3][0] += 2*dz*w;
      A[3][1] += 2*dx*dz*w;
      A[3][2] += 2*dy*dz*w;
      A[3][3] += 2*dz*dz*w;
    }
    //
    matrix_t AI = A.inverse();
    for (int i = 0; i < m_Alpha.size(); ++i) {
      vector_t DX = X[i] - x0;
      value_t  dx = DX[0];
      value_t  dy = DX[1];
      value_t  dz = DX[2];
      m_Alpha[i] = 2*weights[i]*(AI[0][0] + AI[0][1]*dx + AI[0][2]*dy + AI[0][3]*dz);
    }
  }

  int size()
  {
    return m_Alpha.size();
  }

  const std::vector<value_t>& weights() const
  {
    return m_Alpha;
  }

};


} // EDL_NAMESPACE

#include <random>

TEST_CASE("LsqInterpolation3D")
{
  using namespace edl;
  typedef float real;
  typedef MathVector<StaticVector<real,3>> vec_t;
  //
  const int  num_points  = 10;
  const int  num_loops   = 10;
  const int  num_vars    = 5;
  const real range       = 3;
  //
  const std::vector<real> dist_exp = {0, 0.5, 1, 2};
  //
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<real> dist1(-1, 1);
  std::uniform_real_distribution<real> dist2(-range/2, range/2);
  //
  bool all_good = true;
  for (real de : dist_exp) {
    real max_err = 0;
    for (int i_loop = 0; i_loop < num_loops; ++i_loop) {
      //
      // linear test function
      //
      real a[num_vars], b[num_vars], c[num_vars], d[num_vars];
      for (int i = 0; i < num_vars; ++i) {
        a[i] = dist2(gen);
        b[i] = dist2(gen);
        c[i] = dist2(gen);
        d[i] = dist2(gen);
      }
      //
      vec_t g1(dist1(gen), dist1(gen), dist1(gen));
      //vec_t g1(1,0,0);
      g1.normalise();
      vec_t g2 = orthogonalVector(g1);
      //vec_t g2 = vec_t(0,1,0);
      g2.normalise();
      vec_t g3 = g1.cross(g2);
      g3.normalise();
      //
      std::vector<vec_t> x(num_points);
      for (int i = 0; i < num_points; ++i) {
        real x1 = dist1(gen);
        real x2 = dist1(gen);
        real x3 = dist1(gen);
        x[i] = x1*g1 + x2*g2 + x3*g3;
      }
      //
      vec_t x0 = 0.5*g1 + 0.5*g2 + 0.5*g3;
      LsqInterpolation3D<real> lsq;
      lsq.setNodes(x0, x, de);
      //
      for (int i_var = 0; i_var < num_vars; ++i_var) {
        real f0 = a[i_var] + b[i_var]*x0[0] + c[i_var]*x0[1] + d[i_var]*x0[2];
        real fi = 0;
        for (int i = 0; i < num_points; ++i) {
          fi += lsq.weights()[i]*(a[i_var] + b[i_var]*x[i][0] + c[i_var]*x[i][1] + d[i_var]*x[i][2]);
        }
        real err = std::abs(f0 - fi);
        max_err = std::max(max_err, err);
        if (err > 1e-4*range) {
          all_good = false;
        }
      }
    }
  }
  CHECK(all_good);
}

#endif // TGRAPH_H


