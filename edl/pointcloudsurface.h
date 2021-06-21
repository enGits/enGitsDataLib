// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2020 enGits GmbH                                         +
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

#ifndef POINTCLOUDSURFACE_H
#define POINTCLOUDSURFACE_H

#include "edlerror.h"
#include "tlist.h"
#include "mathvector.h"
#include "tsparsetwodimarray.h"
#include "smallsquarematrix.h"
#include "linsolve.h"
#include "geometrytools.h"

#include <thread>
#include <QThread>

namespace EDL_NAMESPACE
{
template <class T> class PointCloudSurface;
}

#include <cstddef>
#include <cstdlib>
#include <typeinfo>
#include <limits>

namespace EDL_NAMESPACE
{

template <class T>  
class PointCloudSurface : public List
{

public: // data types

  typedef MathVector<StaticVector<T,3>> vec_t;
  typedef SmallSquareMatrix<T,3>        mat_t;

  struct weight_t
  {
    T      weight;
    size_t index;
  };

  struct DistanceWeightedAverage
  {
    static void compute(TList<vec_t>& points, TList<vec_t>&, vec_t x, vec_t, TList<T>& weights);
  };

  struct DistanceWeightedLeastSquares
  {
    static void compute(TList<vec_t>& points, TList<vec_t>& normals, vec_t x, vec_t n, TList<T>& weights);
  };


private: // data types

  friend class WeightThread;

  template <class I>
  class WeightThread : public QThread
  {
  public:
    PointCloudSurface<T>* m_PCS;
    std::vector<std::list<weight_t>> m_LocalWeights;
    size_t m_I1, m_I2;
    WeightThread(PointCloudSurface<T>* pcs, size_t i1, size_t i2);
    virtual void run();
  };

  QMutex m_Mutex;

  
protected: // attributes
  
  TList<vec_t>                 m_Points;
  TList<vec_t>                 m_Normals;
  TList<vec_t>                 m_InterpolationPoints;
  TList<vec_t>                 m_InterpolationNormals;
  TSparseTwoDimArray<weight_t> m_Weights;
  bool                         m_SizeDefined = false;
  T                            m_WeightThreshold = 0.01;
  bool                         m_WeightsReady = false;
  vec_t                        m_PlaneOrigin = vec_t(0,0,0);
  vec_t                        m_PlaneNormal = vec_t(1,0,0);
  size_t                       m_MaxStencilSize;
  size_t                       m_MinStencilSize;
  T                            m_AveStencilSize;


public:

  PointCloudSurface();
  ~PointCloudSurface();
  
  template <class C> 
  void setPoints(const C& points);
  
  template <class C> 
  void setNormals(const C& normals);
  void setNormals(const vec_t& normal);
  
  template <class C>
  void setInterpolationPoints(const C& points);

  template <class C>
  void setInterpolationPoints(const C& points, vec_t normal);

  template <class C1, class C2>
  void setInterpolationPoints(const C1& points, const C2& normals);
  
  template <class I> 
  void computeWeights();
  bool weightsReady() { return m_WeightsReady; }

  int numPoints() { return m_Points.size(); }
  int numInterpolationPoints() { return m_InterpolationPoints.size(); }

  template <class C1, class C2>
  void interpolate(const C1& src, C2& dst, int vars_per_point=1);

  void planeFit(T tol=1e-8);
  vec_t planeNormal() { return m_PlaneNormal; }
  vec_t planeOrigin() { return m_PlaneOrigin; }

  int maxInterpolationStencilSize() { return m_MaxStencilSize; }
  int minInterpolationStencilSize() { return m_MinStencilSize; }
  T   averageInterpolationStencilSize() { return m_AveStencilSize; }

  void setWeightsThreshold(T threshold) { m_WeightThreshold = threshold; }

};

template <class T>
PointCloudSurface<T>::PointCloudSurface()
  : List(10, 10), m_InterpolationPoints(0, 10), m_InterpolationNormals(&m_InterpolationPoints), m_Weights(&m_InterpolationPoints)
{
  m_Points.link(this);
  m_Normals.link(this);
  m_Normals.setDefaultValue(vec_t(1,0,0));
}

template <class T>
PointCloudSurface<T>::~PointCloudSurface()
{
  int dummy=0;
}

template <class T>
template <class C>
void PointCloudSurface<T>::setPoints(const C& points)
{
  if (!m_SizeDefined) {
    alloc(points.size());
    m_SizeDefined = true;
  }
  if (points.size() != numActiveEntries()) {
    throw EdlError("size of provided container does not match");
  }
  //
  typename C::const_iterator iter = points.begin();
  size_t i = 0;
  while (iter != points.end()) {
    m_Points[i] = *iter;
    ++iter;
    ++i;
  }
}

template <class T>
template <class C>
void PointCloudSurface<T>::setNormals(const C& normals)
{
  if (!m_SizeDefined) {
    alloc(normals.size());
    m_SizeDefined = true;
  }
  if (normals.size() != numActiveEntries()) {
    throw EdlError("size of provided container does not match");
  }
  //
  typename C::const_iterator iter = normals.begin();
  size_t i = 0;
  while (iter != normals.end()) {
    m_Normals[i] = *iter;
    ++iter;
    ++i;
  }
}

template <class T>
void PointCloudSurface<T>::setNormals(const vec_t& normal)
{
  if (!m_SizeDefined) {
    throw EdlError("please call setPoints first");
  }
  FORALL(i, this->) {
    m_Normals[i] = normal;
  }
}

template <class T>
template <class C>
void PointCloudSurface<T>::setInterpolationPoints(const C& points)
{
  size_t delta = points.size() - m_InterpolationPoints.maxNumEntries();
  if (m_InterpolationPoints.numActiveEntries() > 0) {
    m_InterpolationPoints.delAll();
  }
  m_InterpolationPoints.setDelta(delta);
  m_InterpolationPoints.alloc(points.size());
  //
  typename C::const_iterator iter = points.begin();
  size_t i = 0;
  while (iter != points.end()) {
    m_InterpolationPoints[i] = *iter;
    ++iter;
    ++i;
  }
}

template <class T>
template <class C1, class C2>
void PointCloudSurface<T>::setInterpolationPoints(const C1& points, const C2& normals)
{
  if (points.size() != normals.size()) {
    throw EdlError("points and normals must have the same size");
  }
  //
  size_t delta = points.size() - m_InterpolationPoints.maxNumEntries();
  if (m_InterpolationPoints.numActiveEntries() > 0) {
    m_InterpolationPoints.delAll();
  }
  m_InterpolationPoints.setDelta(delta);
  m_InterpolationPoints.alloc(points.size());
  //
  typename C1::const_iterator iter1 = points.begin();
  typename C2::const_iterator iter2 = normals.begin();
  size_t i = 0;
  while (iter1 != points.end()) {
    m_InterpolationPoints[i] = *iter1;
    m_InterpolationNormals[i] = *iter2;
    ++iter1;
    ++iter2;
    ++i;
  }
}

template <class T>
template <class C>
void PointCloudSurface<T>::setInterpolationPoints(const C& points, vec_t normal)
{
  size_t delta = points.size() - m_InterpolationPoints.maxNumEntries();
  if (m_InterpolationPoints.numActiveEntries() > 0) {
    m_InterpolationPoints.delAll();
  }
  m_InterpolationPoints.setDelta(delta);
  m_InterpolationPoints.alloc(points.size());
  //
  typename C::const_iterator iter = points.begin();
  size_t i = 0;
  while (iter != points.end()) {
    m_InterpolationPoints[i] = *iter;
    m_InterpolationNormals[i] = normal;
    ++iter;
    ++i;
  }
}

template <class T>
template <class I>
void PointCloudSurface<T>::computeWeights()
{
  m_Weights.resetData(m_InterpolationPoints.size(), m_InterpolationPoints.size());
  int num_cores = 1; //std::min(m_InterpolationPoints.size(), std::max(size_t(1), size_t(std::thread::hardware_concurrency())));
  std::vector<WeightThread<I>*> threads(num_cores);
  for (int ith = 0; ith < num_cores; ++ith) {
    size_t i1 = ith*m_InterpolationPoints.size()/num_cores;
    size_t i2 = std::min(m_InterpolationPoints.size(), (ith+1)*m_InterpolationPoints.size()/num_cores);
    threads[ith] = new WeightThread<I>(this, i1, i2);
  }
  for (int ith = 0; ith < num_cores; ++ith) {
    threads[ith]->start();
  }
  for (int ith = 0; ith < num_cores; ++ith) {
    threads[ith]->wait();
    for (size_t i = threads[ith]->m_I1; i < threads[ith]->m_I2; ++i) {
      for (weight_t wgt : threads[ith]->m_LocalWeights[i - threads[ith]->m_I1]) {
        m_Weights.newSubEntry(i) = wgt;
      }
    }
    delete threads[ith];
  }
  //
  m_MinStencilSize = m_Weights.size();
  m_MaxStencilSize = 0;
  m_AveStencilSize = 0;
  int N = 0;
  FORALL(i, m_Weights.) {
    ++N;
    m_MinStencilSize = std::min(m_MinStencilSize, m_Weights.count(i));
    m_MaxStencilSize = std::max(m_MaxStencilSize, m_Weights.count(i));
    m_AveStencilSize += T(m_Weights.count(i));
  }
  m_AveStencilSize /= N;
  //
  m_WeightsReady = true;
}

template <class T>
template <class C1, class C2>
void PointCloudSurface<T>::interpolate(const C1& src, C2& dst, int vars_per_point)
{
  T _zero = 0;
  typename C2::value_type zero(_zero);

//  FORALL(i, m_InterpolationPoints.) {
//    dst[i] = zero;
//    for (int j = 0; j < m_Weights.count(i); ++j) {
//      size_t idx = m_Weights.at(i,j).index;
//      dst[i] += m_Weights.at(i,j).weight * src[idx];
//    }
//  }

  FORALL (i_point, m_InterpolationPoints.) {
    for (int i_var = 0; i_var < vars_per_point; ++i_var) {
      size_t dst_idx = i_point*vars_per_point + i_var;
      dst[dst_idx] = zero;
      for (int i_contrib = 0; i_contrib < m_Weights.count(i_point); ++i_contrib) {
        size_t src_idx = m_Weights.at(i_point, i_contrib).index*vars_per_point + i_var;
        dst[dst_idx] += m_Weights.at(i_point, i_contrib).weight * src[src_idx];
      }
    }
  }
}

template <class T>
void PointCloudSurface<T>::planeFit(T tol)
{
  SmallSquareMatrix<T,4> A;
  A.initAll(0);
  vec_t sum(0,0,0);
  vec_t sum_sq(0,0,0);
  FORALL(i, m_Points.) {
    //
    A[0][0] += 1;
    A[0][1] -= m_Points[i][0];
    A[0][2] -= m_Points[i][1];
    A[0][3] -= m_Points[i][2];
    //
    A[1][0] -= m_Points[i][0];
    A[1][1] += m_Points[i][0]*m_Points[i][0];
    A[1][2] += m_Points[i][0]*m_Points[i][1];
    A[1][3] += m_Points[i][0]*m_Points[i][2];
    //
    A[2][0] -= m_Points[i][1];
    A[2][1] += m_Points[i][1]*m_Points[i][0];
    A[2][2] += m_Points[i][1]*m_Points[i][1];
    A[2][3] += m_Points[i][1]*m_Points[i][2];
    //
    A[3][0] -= m_Points[i][2];
    A[3][1] += m_Points[i][2]*m_Points[i][0];
    A[3][2] += m_Points[i][2]*m_Points[i][1];
    A[3][3] += m_Points[i][2]*m_Points[i][2];
    //
    for (int j = 0; j < 3; ++j) {
      sum[j]    += m_Points[i][j];
      sum_sq[j] += m_Points[i][j]*m_Points[i][j];
    }
  }

  vec_t x_mean, x_std, n;
  for (int i = 0; i < 3; ++i) {
    x_mean[i] = sum[i]/m_Points.size();
    x_std[i]  = sqrt(sum_sq[i]/m_Points.size() - x_mean[i]*x_mean[i]);
    n[i]      = 1.0/std::max(T(0.01), x_std[i]);
  }
  n.normalise();

  MathVector<StaticVector<T,4>> rhs(0,0,0,0), result(0, n[0], n[1], n[2]);
  try {
    iterLinsolve(A, rhs, result, tol);
  } catch (LinSolveError e) {
    std::cerr << "cannot fit plane (det(A)=" << e.det << ")" << std::endl;
    throw e;
  }
  n[0] = result[1];
  n[1] = result[2];
  n[2] = result[3];
  n.normalise();
  m_PlaneOrigin = result[0]*n;
  m_PlaneNormal = n;
}

template <class T>
template <class I>
PointCloudSurface<T>::WeightThread<I>::WeightThread(PointCloudSurface<T>* pcs, size_t i1, size_t i2)
{
  m_PCS = pcs;
  m_I1  = i1;
  m_I2  = i2;
}

template <class T>
template <class I>
void PointCloudSurface<T>::WeightThread<I>::run()
{
  m_PCS->m_Mutex.lock();
  TList<T>* global_weights = new TList<T>(m_PCS);
  m_PCS->m_Mutex.unlock();
  m_LocalWeights.resize(m_I2 - m_I1, std::list<weight_t>());
  for (size_t i = m_I1; i < m_I2; ++i) {
    I::compute(m_PCS->m_Points, m_PCS->m_Normals, m_PCS->m_InterpolationPoints[i], m_PCS->m_InterpolationNormals[i], *global_weights);
    T total_weight = 0;
    T max_weight = 0;
    FORALL(j, global_weights->) {
      max_weight = std::max(max_weight, global_weights->at(j));
      if (global_weights->at(j) > m_PCS->m_WeightThreshold) {
        total_weight += global_weights->at(j);
        weight_t wgt;
        wgt.index = j;
        wgt.weight = global_weights->at(j);
        m_LocalWeights[i - m_I1].push_back(wgt);
      }
    }
    for (weight_t& wgt : m_LocalWeights[i - m_I1]) {
      wgt.weight /= total_weight;
    }
  }
  m_PCS->m_Mutex.lock();
  delete global_weights;
  m_PCS->m_Mutex.unlock();
}




template <class T>
void PointCloudSurface<T>::DistanceWeightedAverage::compute(TList<vec_t>& points, TList<vec_t>&, vec_t x, vec_t, TList<T>& weights)
{
  T min_d = std::numeric_limits<T>::max();
  std::list<T> min_ds;
  bool sorted = false;
  FORALL(i, weights.) {    
    T d = (points[i] - x).abs();
    if (min_ds.size() < 4) {
      min_ds.push_back(d);
    } else {
      if (!sorted) {
        min_ds.sort();
        sorted = true;
      }
      if (d < min_ds.back()) {
        min_ds.pop_back();
        min_ds.push_front(d);
        min_ds.sort();
      }
    }
    min_d = std::min(min_d, d);
    weights[i] = d;
  }
  //
  min_ds.sort();
  min_d = std::max(min_d, std::numeric_limits<T>::min());
  T total_weight = 0;
  //
  FORALL(i, points.) {
    if (weights[i] <= min_ds.back()) {
      weights[i] = 1.0/std::max(min_d, weights[i]);
      total_weight += weights[i];
    } else {
      weights[i] = 0;
    }
  }
  //
  FORALL(i, points.) {
    weights[i] /= total_weight;
  }
}
  
template <class T>
void PointCloudSurface<T>::DistanceWeightedLeastSquares::compute(TList<vec_t>& points, TList<vec_t>& normals, vec_t x, vec_t n, TList<T>& weights)
{
  vec_t g1 = n;
  vec_t g2 = orthogonalVector(n);
  vec_t g3 = g1.cross(g2);
  g1.normalise();
  g2.normalise();
  g3.normalise();
  mat_t G;
  G.setColumn(0, g1);
  G.setColumn(1, g2);
  G.setColumn(2, g3);
  mat_t GI = G.inverse();
  //
  vec_t x2(1,0,0);
  //
  T min_dist = std::numeric_limits<T>::max();
  FORALL(i, weights.) {
    vec_t dx = points[i] - x;
    T dxn = dx*g1;
    dx -= dxn*g1;
    //weights[i] = (points[i] - x).abs();
    weights[i] = dx.abs();
    min_dist = std::min(min_dist, weights[i]);
  }
  //
  min_dist = std::max(min_dist, std::numeric_limits<T>::min());
  //
  mat_t A;
  A.initAll(0);
  FORALL(i, points.) {
    vec_t xi2 = points[i] - x;
    xi2 = GI*xi2;
    xi2[0] = 1;
    weights[i] = 1.0/std::max(min_dist, weights[i]);
    A[0] += weights[i]*weights[i] *        xi2;
    A[1] += weights[i]*weights[i] * xi2[1]*xi2;
    A[2] += weights[i]*weights[i] * xi2[2]*xi2;
  }
  mat_t AI = A.inverse();
  T total_weight = 0;
  FORALL(i, points.) {
    vec_t xi2 = points[i] - x;
    xi2 = GI*xi2;
    xi2[0] = 1;
    weights[i] = weights[i]*weights[i] * (AI*xi2)*x2;
    total_weight += weights[i];
  }
  FORALL(i, points.) {
    weights[i] /= total_weight;
  }
}
  
  
} // namespace

#endif // POINTCLOUDSURFACE_H

