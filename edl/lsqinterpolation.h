// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef LSQINTERPOLATION_H
#define LSQINTERPOLATION_H

#include "edl/edl.h"
#include "edl/mathvector.h"
#include "edl/smallsquarematrix.h"

#include <iostream>
#include <sys/types.h>
#include <vector>

namespace EDL_NAMESPACE
{
template <typename T> class LsqInterpolation3D;
template <typename T> class LsqInterpolationCoPlanar;
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


protected: // methods

  void computeWeights()
  {
    auto N = m_CollectedNodes.size();
    m_Alpha.resize(N);
    matrix_t A;
    A.initAll(0);
    for (int i = 0; i < N; ++i) {
      geovec_t DX = m_CollectedNodes[i] - m_X0;
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
      geovec_t DX = m_CollectedNodes[i] - m_X0;
      value_t  dx = DX[0];
      value_t  dy = DX[1];
      value_t  dz = DX[2];
      m_Alpha[i] = 2*m_CollectedWeights[i]*(AI[0][0] + AI[0][1]*dx + AI[0][2]*dy + AI[0][3]*dz);
    }
    m_Valid = true;
  }


public: // methods

  LsqInterpolation3D(geovec_t x0 = geovec_t(0,0,0)) : m_X0(x0) {}

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
      value_t d = std::max((x - m_X0).abs(), min_dist);
      value_t b = d/min_dist;
      value_t w = std::pow(b, distance_exponent);
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

};

template <typename VALUE_TYPE>
class LsqInterpolationCoPlanar : public LsqInterpolation3D<VALUE_TYPE>
{

protected: // data types

  using typename LsqInterpolation3D<VALUE_TYPE>::value_t;
  using typename LsqInterpolation3D<VALUE_TYPE>::geovec_t;

  typedef MathVector<StaticVector<value_t,3>> vector_t;
  typedef SmallSquareMatrix<value_t,3>        matrix_t;


protected: // methods

  void projectToPlane()
  {
    geovec_t x0 = origin();
    geovec_t n  = normal();
    geovec_t g3 = n;
    geovec_t g1 = orthogonalVector(g3);
    geovec_t g2 = g3.cross(g1);
    //
    // base matrix
    //
    matrix_t B;
    B.initAll(0);
    B.setColumn(0, g1);
    B.setColumn(1, g2);
    B.setColumn(2, g3);
    matrix_t BI = B.inverse();
    //
    // project nodes
    //
    for (int i = 0; i < this->m_CollectedNodes.size(); ++i) {
      geovec_t dx = (this->m_CollectedNodes[i] - x0);
      this->m_CollectedNodes[i] = BI*dx;
    }
    //
    // project origin
    //
    geovec_t dx = this->m_X0 - x0;
    this->m_X0 = BI*dx;
  }

  void computeWeights()
  {
    projectToPlane();
    auto N = this->m_CollectedNodes.size();
    this->m_Alpha.resize(N);
    matrix_t A;
    A.initAll(0);
    for (int i = 0; i < N; ++i) {
      vector_t DX = this->m_CollectedNodes[i] - this->m_X0;
      value_t  dx = DX[0];
      value_t  dy = DX[1];
      value_t  dz = DX[2];
      //
      auto w = this->m_CollectedWeights[i];
      //
      A[0][0] += 2*w;
      A[0][1] += 2*dx*w;
      A[0][2] += 2*dy*w;
      A[1][0] += 2*dx*w;
      A[1][1] += 2*dx*dx*w;
      A[1][2] += 2*dx*dy*w;
      A[2][0] += 2*dy*w;
      A[2][1] += 2*dx*dy*w;
      A[2][2] += 2*dy*dy*w;
    }
    //
    matrix_t AI = A.inverse();
    for (int i = 0; i < this->m_Alpha.size(); ++i) {
      vector_t DX = this->m_CollectedNodes[i] - this->m_X0;
      value_t  dx = DX[0];
      value_t  dy = DX[1];
      value_t  dz = DX[2];
      this->m_Alpha[i] = 2*this->m_CollectedWeights[i]*(AI[0][0] + AI[0][1]*dx + AI[0][2]*dy);
    }
    this->m_Valid = true;
  }


public:

  LsqInterpolationCoPlanar(geovec_t x0 = geovec_t(0,0,0)) : LsqInterpolation3D<VALUE_TYPE>(x0) {}

  const std::vector<value_t>& weights()
  {
    if (!this->m_Valid) {
      computeWeights();
    }
    return this->m_Alpha;
  }  

  geovec_t origin()
  {
    geovec_t x(0,0,0);
    for (auto _x : this->m_CollectedNodes) {
      x += _x;
    }
    return value_t(1)/this->m_CollectedNodes.size() * x;
  }

  geovec_t normal()
  {
    geovec_t x0 = origin();
    geovec_t n(0,0,0);
    for (int i = 0; i < this->m_CollectedNodes.size(); ++i) {
      geovec_t dxi = this->m_CollectedNodes[i] - x0;
      for (int j = i+1; j < this->m_CollectedNodes.size(); ++j) {
        geovec_t dxj = this->m_CollectedNodes[j] - x0;
        geovec_t dn  = dxi.cross(dxj);
        if (dn*n < 0) {
          dn *= -1;
        }
        n += dn;
      }
    }
    return n.normalised();
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


TEST_CASE("LsqInterpolationCoPlanar (one interpolation only)")
{
  using namespace edl;
  typedef float real;
  typedef MathVector<StaticVector<real,3>> vec_t;
  //
  vec_t g1 = vec_t(1,0,0);
  g1.normalise();
  vec_t g2 = vec_t(0,0,1);//orthogonalVector(g1);
  g2.normalise();
  vec_t g3 = g1.cross(g2);
  std::vector<vec_t> X;
  X.push_back(vec_t(0,0,0));
  X.push_back(1*g1 + 0*g2 + 0*g3);
  X.push_back(0*g1 + 1*g2 + 0*g3);
  X.push_back(1*g1 + 1*g2 + 0*g3);
  vec_t x0 = 0.5*g1 + 0.5*g2;
  //
  real a = 1;
  real b = 2;
  real c = 3;
  real d = 4;
  real f1 = a + b*x0[0] + c*x0[1] + d*x0[2];
  //
  LsqInterpolationCoPlanar<real> lsq(x0);
  lsq.setNodes(X);
  auto n = lsq.normal();
  if (n*g3 < 0) {
    n *= -1;
  }
  CHECK(n[0] == doctest::Approx(g3[0]));
  CHECK(n[1] == doctest::Approx(g3[1]));
  CHECK(n[2] == doctest::Approx(g3[2]));
  auto w = lsq.weights();
  real f2 = 0.0;
  for (int i = 0; i < w.size(); ++i) {
    f2 += w[i]*(a + b*X[i][0] + c*X[i][1] + d*X[i][2]);
  }
  //
  CHECK(f1 == doctest::Approx(f2));
}


TEST_CASE("LsqInterpolationCoPlanar (setNodes)")
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
      real a[num_vars], b[num_vars], c[num_vars];
      for (int i = 0; i < num_vars; ++i) {
        a[i] = dist2(gen);
        b[i] = dist2(gen);
        c[i] = dist2(gen);
      }
      //
      vec_t g1(dist1(gen), dist1(gen), dist1(gen));
      //vec_t g1(1,0,0);
      g1.normalise();
      vec_t g2 = orthogonalVector(g1);
      //vec_t g2 = vec_t(0,1,0);
      g2.normalise();
      //
      std::vector<vec_t> x;
      x.reserve(num_points);
      for (int i = 0; i < num_points - 1; ++i) {
        real x1 = dist1(gen);
        real x2 = dist1(gen);
        x.push_back(x1*g1 + x2*g2);
      }
      //
      vec_t x0 = 0.5*g1 + 0.5*g2;
      x.push_back(x0);
      LsqInterpolationCoPlanar<real> lsq(x0);
      lsq.setNodes(x, de);
      //
      for (int i_var = 0; i_var < num_vars; ++i_var) {
        real f0 = a[i_var] + b[i_var]*x0[0] + c[i_var]*x0[1];
        real fi = 0;
        for (int i = 0; i < num_points; ++i) {
          fi += lsq.weights()[i]*(a[i_var] + b[i_var]*x[i][0] + c[i_var]*x[i][1]);
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


#endif // LSQINTERPOLATION_H


