// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef POINTCLOUDSURFACE_H
#define POINTCLOUDSURFACE_H

#include "edl/edl.h"
#include "edl/edlerror.h"
#include "edl/mathvector.h"
#include "edl/lsqinterpolation.h"
#include "edl/geometrytools.h"
#include "edl/smallsquarematrix.h"

#include <array>
#include <iostream>
#include <vector>

namespace EDL_NAMESPACE
{
template <class T, uint8_t NUM_WEIGHTS> class PointCloudSurface;
}

#include <limits>

namespace EDL_NAMESPACE
{

template <class T, uint8_t NUM_WEIGHTS>  
class PointCloudSurface
{

public: // data types

  typedef MathVector<StaticVector<T,3>> vec_t;
  typedef SmallSquareMatrix<T,3>        mat_t;
  // typedef LsqInterpolation3D<T>         lsq_t;

  struct weights_t
  {
    T      weight[NUM_WEIGHTS];
    size_t index [NUM_WEIGHTS];
  };

private: // data types

protected: // attributes
  
  std::vector<vec_t>     m_Points;
  std::vector<vec_t>     m_InterpolationPoints;
  std::vector<weights_t> m_Weights;
  T                      m_DistExp;
  vec_t                  m_G1;
  vec_t                  m_G2;
  vec_t                  m_G3;
  vec_t                  m_X0;


public:

  PointCloudSurface(vec_t n, T dist_exp = 0) : m_G3(n), m_DistExp(dist_exp) 
  {
    m_G3.normalise();
    m_G2 = orthogonalVector(m_G3);
    m_G1 = m_G2.cross(m_G3);
  }

  template <class C> 
  void setPoints(const C& points)
  {
    m_Points.clear();
    m_Points.reserve(points.size());
    for (const vec_t& p : points) {
      m_Points.push_back(p);
    }
  }
  
  template <class C>
  void setInterpolationPoints(const C& points)
  {
    m_InterpolationPoints.clear();
    m_InterpolationPoints.reserve(points.size());
    for (const vec_t& p : points) {
      m_InterpolationPoints.push_back(p);
    }
  }
  
  void computeWeights()
  {
    using std::vector;
    m_Weights.clear();
    m_Weights.resize(m_InterpolationPoints.size());
    //
    m_X0 = vec_t(0,0,0);
    for (vec_t x : m_Points) {
      m_X0 += x;
    }
    m_X0 *= T(1)/m_Points.size();
    //
    mat_t M;
    M.setColumn(0, m_G1);
    M.setColumn(1, m_G2);
    M.setColumn(2, m_G3);
    M = M.inverse();
    //
    for (size_t i = 0; i < m_InterpolationPoints.size(); ++i) {
      vec_t x0 = m_InterpolationPoints[i];
      //
      // find the NUM_WEIGHTS closest points
      //
      vector<std::pair<T,size_t>> dists;
      dists.reserve(m_Points.size());
      for (size_t j = 0; j < m_Points.size(); ++j) {
        vec_t xj = m_Points[j];
        T dist = (x0-xj).abs();
        dists.push_back(std::make_pair(dist, j));
      }
      std::sort(dists.begin(), dists.end());
      //
      LsqInterpolationCoPlanar<T> lsq;
      //
      size_t num_weights = std::min(size_t(NUM_WEIGHTS), m_Points.size());
      std::vector<vec_t> X(num_weights);
      for (int j = 0; j < num_weights; ++j) {
        m_Weights[i].index[j] = dists[j].second;
        vec_t x3d = m_Points[dists[j].second] - x0;
        vec_t x2d = M*x3d;
        X[j] = vec_t(x2d[0], x2d[1], T(0));
      }
      lsq.setNodes(X, m_DistExp);
      for (int j = 0; j < num_weights; ++j) {
        m_Weights[i].weight[j] = lsq.weights()[j];
      }
      for (int j = num_weights; j < NUM_WEIGHTS; ++j) {
        m_Weights[i].weight[j] = T(0);
      }
    }
  }

  int numPoints() { return m_Points.size(); }
  int numInterpolationPoints() { return m_InterpolationPoints.size(); }

  template <class C1, class C2>
  void interpolate(const C1& src, C2& dst, int vars_per_point=1)
  {
    for (int i = 0; i < m_InterpolationPoints.size(); ++i) {
      for (int j = 0; j < vars_per_point; ++j) {
        dst[i*vars_per_point + j] = T(0);
        for (int k = 0; k < NUM_WEIGHTS; ++k) {
          size_t idx = m_Weights[i].index[k];
          real   w   = m_Weights[i].weight[k];
          dst[i*vars_per_point + j] += w*src[idx*vars_per_point + j];
        }
      }
    }
  }

};

} // namespace

#include <random>

TEST_CASE("PointCloudSurface")
{
  using namespace edl;
  typedef float real;
  typedef MathVector<StaticVector<real,3>> vec_t;
  typedef SmallSquareMatrix<real,3>        mat_t;
  //
  const int  size_i      = 100;
  const int  size_j      = 100;
  const int  num_ipoints = 50;
  const int  num_loops   = 10;
  const int  num_vars    = 5;
  const real range       = 3;
  //
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<real> dist1(0.1, 1.0);
  std::uniform_real_distribution<real> dist2(0, 1);
  std::uniform_real_distribution<real> dist3(-3, 3);
  //
  bool all_good = true;
  real err_max  = 0;
  for (int i_loop = 0; i_loop < num_loops; ++i_loop) {
    //
    vec_t g1(dist1(gen), dist1(gen), dist1(gen));
    g1.normalise();
    vec_t g2 = orthogonalVector(g1);
    g2.normalise();
    vec_t g3 = g1.cross(g2);
    g3.normalise();
    //
    // linear test function
    //
    mat_t T;
    T[0] = g1;
    T[1] = g2;
    T[2] = g3;
    T = T.transp();
    real a[num_vars], b[num_vars], c[num_vars], d[num_vars];
    for (int i = 0; i < num_vars; ++i) {
      vec_t _C(0, dist3(gen), dist3(gen));
      vec_t C = T*_C;
      a[i] = dist3(gen);
      b[i] = C[0];
      c[i] = C[1];
      d[i] = C[2];
    }
    //
    std::vector<vec_t> points(size_i*size_j);
    for (int i = 0; i < size_i; ++i) {
      for (int j = 0; j < size_j; ++j) {
        real c2 = i*1.0/(size_i-1);
        real c3 = j*1.0/(size_j-1);
        points[i*size_j + j] = c2*g2 + c3*g3;
      }
    }
    //
    std::vector<vec_t> ipoints(num_ipoints);
    for (int i = 0; i < num_ipoints; ++i) {
      real c1 = 0.05*dist2(gen);
      real c2 = dist2(gen);
      real c3 = dist2(gen);
      ipoints[i] = c1*g1 + c2*g2 + c3*g3;
    }
    //
    std::vector<real> vars(num_vars*points.size());
    for (int i = 0; i < points.size(); ++i) {
      for (int j = 0; j < num_vars; ++j) {
        vars[i*num_vars + j] = a[j] + b[j]*points[i][0] + c[j]*points[i][1] + d[j]*points[i][2];
      }
    }
    //
    std::vector<real> ivars(num_vars*ipoints.size());
    for (int i = 0; i < ipoints.size(); ++i) {
      for (int j = 0; j < num_vars; ++j) {
        ivars[i*num_vars + j] = a[j] + b[j]*ipoints[i][0] + c[j]*ipoints[i][1] + d[j]*ipoints[i][2];
      }
    }
    //
    PointCloudSurface<real,10> pcs(g1, 0);
    pcs.setPoints(points);
    pcs.setInterpolationPoints(ipoints);
    pcs.computeWeights();
    std::vector<real> result(num_vars*ipoints.size());
    pcs.interpolate(vars, result, num_vars);
    //
    for (int i = 0; i < ipoints.size(); ++i) {
      for (int j = 0; j < num_vars; ++j) {
        real f0 = ivars[i*num_vars + j];
        real f1 = result[i*num_vars + j];
        real err = std::abs(f0-f1);
        err_max = std::max(err_max, err);
        if (err > 1e-4*range) {
          all_good = false;
        }
      }
    }
  }
  CHECK(all_good);
}

#endif // POINTCLOUDSURFACE_H

