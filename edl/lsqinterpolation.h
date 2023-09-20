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

#include <sys/types.h>
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


public: // data types

  enum class mode_t
  {
    XYZ,
    XY,
    XZ,
    YZ    
  };


protected: // attributes

  geovec_t              m_X0;
  std::vector<value_t>  m_Alpha;
  std::vector<value_t>  m_CollectedWeights;
  std::vector<geovec_t> m_CollectedNodes;
  bool                  m_Valid{false};
  mode_t                m_Mode{mode_t::XYZ};


protected: // methods

  void computeWeights()
  {
    auto N = m_CollectedNodes.size();
    m_Alpha.resize(N);
    matrix_t A;
    A.initAll(0);
    for (int i = 0; i < N; ++i) {
      vector_t DX = m_CollectedNodes[i] - m_X0;
      value_t  dx = DX[0];
      value_t  dy = DX[1];
      value_t  dz = DX[2];
      //
      auto w = m_CollectedWeights[i];
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
      vector_t DX = m_CollectedNodes[i] - m_X0;
      value_t  dx = DX[0];
      value_t  dy = DX[1];
      value_t  dz = DX[2];
      m_Alpha[i] = 2*m_CollectedWeights[i]*(AI[0][0] + AI[0][1]*dx + AI[0][2]*dy + AI[0][3]*dz);
    }
    m_Valid = true;
  }


public: // methods

  LsqInterpolation3D(geovec_t x0 = geovec_t(0,0,0)) : m_X0(x0) {}

  void setMode(mode_t mode) { m_Mode = mode; }

  template<typename C>
  void setNodes(C& X, value_t distance_exponent=0.0)
  {
    value_t min_dist = std::numeric_limits<value_t>::max();
    value_t max_dist = 0;
    for (auto x : X) {
      min_dist = std::min(min_dist, (x - m_X0).abs());
      max_dist = std::max(max_dist, (x - m_X0).abs());
    }
    min_dist = std::max(min_dist, max_dist/100);
    m_CollectedNodes.clear();
    m_CollectedNodes.reserve(X.size());
    m_CollectedWeights.clear();
    m_CollectedWeights.reserve(X.size());
    for (auto x : X) {
      value_t  w  = std::pow((x - m_X0).abs()/min_dist, distance_exponent);
      addNode(x, w);
    }
  }

  void addNode(geovec_t x, value_t weight)
  {
    m_CollectedNodes.push_back(x);
    m_CollectedWeights.push_back(weight);
    m_Valid = false;
  }

  int size()
  {
    return m_Alpha.size();
  }

  const std::vector<value_t>& weights()
  {
    if (!m_Valid) {
      computeWeights();
    }
    return m_Alpha;
  }

  /**
    * @brief compute a pure interpolation without any extrapolation if possible
    * @return a vector of interpolation weights (>0) if possible
    */
  std::vector<value_t> pureInterpolationWeights()
  {
    using namespace std;
    //
    LsqInterpolation3D<value_t> lsq(m_X0);
    vector<geovec_t> X = m_CollectedNodes;
    vector<value_t>  w;
    lsq.setMode(m_Mode);    
    //
    bool all_positive = false;
    while (!all_positive) {
      lsq.setNodes(X, 0);
      w = lsq.weights();
      vector<pair<int, value_t>> wi(w.size());
      // sort weights
      for (int i = 0; i < w.size(); ++i) {
        wi[i] = make_pair(i, w[i]);
      }
      sort(wi.begin(), wi.end(), [](const pair<int, value_t> &a, const pair<int, value_t> &b) { return a.second < b.second; });
      // check if all weights are positive
      all_positive = wi[0].second > 0;
      if (!all_positive) {
        vector<geovec_t> X_new;
        for (int i = 1; i < wi.size(); ++i) {
          X_new.push_back(m_CollectedNodes[wi[i].first]);
        }
        X = X_new;
      }
    }
    return w;
  }

};

} // EDL_NAMESPACE

// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

#include <random>

TEST_CASE("LsqInterpolation3D (setNodes)")
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
      std::vector<vec_t> x;
      x.reserve(num_points);
      for (int i = 0; i < num_points - 1; ++i) {
        real x1 = dist1(gen);
        real x2 = dist1(gen);
        real x3 = dist1(gen);
        x.push_back(x1*g1 + x2*g2 + x3*g3);
      }
      //
      vec_t x0 = 0.5*g1 + 0.5*g2 + 0.5*g3;
      x.push_back(x0);
      LsqInterpolation3D<real> lsq(x0);
      lsq.setNodes(x, de);
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

TEST_CASE("LsqInterpolation3D (multiple addNode)")
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
      LsqInterpolation3D<real> lsq(x0);
      for (auto _x : x) {
        vec_t Dx = _x - x0;
        real  w  = 1.0/std::max(real(1e-4), Dx.abs());
        lsq.addNode(_x, w);
      }
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

TEST_CASE("LsqInterpolation3D (pureInterpolation)")
{
  using namespace edl;
  using namespace std;
  typedef float real;
  typedef MathVector<StaticVector<real,3>> vec_t;
  //
  vector<real> coeff = {1, 2, 3, 4};
  vector<vec_t> X;
  X.push_back(vec_t(0,0,0));
  X.push_back(vec_t(1,0,0));
  X.push_back(vec_t(2,0,0));
  X.push_back(vec_t(0,1,0));
  X.push_back(vec_t(1,1,0));
  X.push_back(vec_t(2,1,0));
  X.push_back(vec_t(0,0,1));
  X.push_back(vec_t(1,0,1));
  X.push_back(vec_t(2,0,1));
  X.push_back(vec_t(0,1,1));
  X.push_back(vec_t(1,1,1));
  X.push_back(vec_t(2,1,1));
  //
  vector<real> F;
  for (auto x : X) {
    real f = coeff[0] + coeff[1]*x[0] + coeff[2]*x[1] + coeff[3]*x[2];
    F.push_back(f);
  }
  //
  vec_t x0(1.5,0.5,0.5);
  LsqInterpolation3D<real> lsq(x0);
  lsq.setNodes(X, 0);
  auto w = lsq.pureInterpolationWeights();
  //auto w = lsq.weights();
  //
  real f1 = coeff[0] + coeff[1]*x0[0] + coeff[2]*x0[1] + coeff[3]*x0[2];
  real f2 = 0;
  for (int i = 0; i < X.size(); ++i) {
    f2 += w[i]*F[i];
  }
  CHECK(f1 == doctest::Approx(f2));
  //
  // CHECK(w[0]  == doctest::Approx(0.000));
  // CHECK(w[1]  == doctest::Approx(0.125));
  // CHECK(w[2]  == doctest::Approx(0.125));
  // CHECK(w[3]  == doctest::Approx(0.000));
  // CHECK(w[4]  == doctest::Approx(0.125));
  // CHECK(w[5]  == doctest::Approx(0.125));
  // CHECK(w[6]  == doctest::Approx(0.000));
  // CHECK(w[7]  == doctest::Approx(0.125));
  // CHECK(w[8]  == doctest::Approx(0.125));
  // CHECK(w[9]  == doctest::Approx(0.000));
  // CHECK(w[10] == doctest::Approx(0.125));
  // CHECK(w[11] == doctest::Approx(0.125));
}

#endif // LSQINTERPOLATION_H


