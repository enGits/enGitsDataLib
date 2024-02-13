// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2023 enGits GmbH                                         +
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
#ifndef TLSQGRAD_H
#define TLSQGRAD_H

#ifndef CUDA_DH
#define CUDA_DH
#endif

#include "edl/edl.h"


namespace EDL_NAMESPACE
{
template <typename T, uint8_t NUM_VARS_TLSQ> class TLSqGrad3D;
}

namespace EDL_NAMESPACE
{

template <typename T, uint8_t NUM_VARS_TLSQ>
class TLSqGrad3D
{

  T m_Axx{0};
  T m_Axy{0};
  T m_Axz{0};
  T m_Ayx{0};
  T m_Ayy{0};
  T m_Ayz{0};
  T m_Azx{0};
  T m_Azy{0};
  T m_Azz{0};
  T m_Bx[NUM_VARS_TLSQ]{0};
  T m_By[NUM_VARS_TLSQ]{0};
  T m_Bz[NUM_VARS_TLSQ]{0};
  T m_X0{0};
  T m_Y0{0};
  T m_Z0{0};
  T m_F0[NUM_VARS_TLSQ]{0};

  uint16_t m_NumPoints{0};

public:

  CUDA_DH TLSqGrad3D(T x, T y, T z, const T* f) : m_X0(x), m_Y0(y), m_Z0(z)
  {
    for (int i = 0; i < NUM_VARS_TLSQ; ++i) {
      m_F0[i] = f[i];
    }
  }

  CUDA_DH TLSqGrad3D(T x, T y, T z, T f) : m_X0(x), m_Y0(y), m_Z0(z) 
  {
    m_F0[0] = f;
  }

  CUDA_DH void addPoint(T x, T y, T z, const T* f, T w = 1.0)
  {
    T Dx = x - m_X0;
    T Dy = y - m_Y0;
    T Dz = z - m_Z0;
    //
    m_Axx += Dx*Dx*w;
    m_Axy += Dx*Dy*w;
    m_Axz += Dx*Dz*w;
    m_Ayx += Dx*Dy*w;
    m_Ayy += Dy*Dy*w;
    m_Ayz += Dy*Dz*w;
    m_Azx += Dx*Dz*w;
    m_Azy += Dy*Dz*w;
    m_Azz += Dz*Dz*w;
    //
    for (int i = 0; i < NUM_VARS_TLSQ; ++i) {
      T Df = f[i] - m_F0[i];
      m_Bx[i] += Df*Dx*w;
      m_By[i] += Df*Dy*w;
      m_Bz[i] += Df*Dz*w;
    }
    //
    ++m_NumPoints;
  }

  CUDA_DH void addPoint(T x, T y, T z, T f, T w = 1.0)
  {
    addPoint(x, y, z, &f, w);
  }

  CUDA_DH void computeGrad(T (& grad)[NUM_VARS_TLSQ][3])
  {
    T inv_det = T(1) / (m_Axx * m_Ayy * m_Azz - m_Axx * m_Ayz * m_Azy -
                        m_Axy * m_Ayx * m_Azz + m_Axy * m_Ayz * m_Azx +
                        m_Axz * m_Ayx * m_Azy - m_Axz * m_Ayy * m_Azx);
    //
    T AI_xx =  m_Ayy*m_Azz - m_Ayz*m_Azy;
    T AI_xy = -m_Axy*m_Azz + m_Axz*m_Azy;
    T AI_xz =  m_Axy*m_Ayz - m_Axz*m_Ayy;
    T AI_yx = -m_Ayx*m_Azz + m_Ayz*m_Azx;
    T AI_yy =  m_Axx*m_Azz - m_Axz*m_Azx;
    T AI_yz = -m_Axx*m_Ayz + m_Axz*m_Ayx;
    T AI_zx =  m_Ayx*m_Azy - m_Ayy*m_Azx;
    T AI_zy = -m_Axx*m_Azy + m_Axy*m_Azx;
    T AI_zz =  m_Axx*m_Ayy - m_Axy*m_Ayx;
    //
    for (int i = 0; i < NUM_VARS_TLSQ; ++i) {
      grad[i][0] = inv_det*(AI_xx*m_Bx[i] + AI_xy*m_By[i] + AI_xz*m_Bz[i]);
      grad[i][1] = inv_det*(AI_yx*m_Bx[i] + AI_yy*m_By[i] + AI_yz*m_Bz[i]);
      grad[i][2] = inv_det*(AI_zx*m_Bx[i] + AI_zy*m_By[i] + AI_zz*m_Bz[i]);
    }
  }

  CUDA_DH void computeGrad(T* grad)
  {
    T inv_det = T(1) / (m_Axx * m_Ayy * m_Azz - m_Axx * m_Ayz * m_Azy -
                        m_Axy * m_Ayx * m_Azz + m_Axy * m_Ayz * m_Azx +
                        m_Axz * m_Ayx * m_Azy - m_Axz * m_Ayy * m_Azx);
    //
    T AI_xx =  m_Ayy*m_Azz - m_Ayz*m_Azy;
    T AI_xy = -m_Axy*m_Azz + m_Axz*m_Azy;
    T AI_xz =  m_Axy*m_Ayz - m_Axz*m_Ayy;
    T AI_yx = -m_Ayx*m_Azz + m_Ayz*m_Azx;
    T AI_yy =  m_Axx*m_Azz - m_Axz*m_Azx;
    T AI_yz = -m_Axx*m_Ayz + m_Axz*m_Ayx;
    T AI_zx =  m_Ayx*m_Azy - m_Ayy*m_Azx;
    T AI_zy = -m_Axx*m_Azy + m_Axy*m_Azx;
    T AI_zz =  m_Axx*m_Ayy - m_Axy*m_Ayx;
    //
    for (int i = 0; i < NUM_VARS_TLSQ; ++i) {
      grad[0] = inv_det*(AI_xx*m_Bx[i] + AI_xy*m_By[i] + AI_xz*m_Bz[i]);
      grad[1] = inv_det*(AI_yx*m_Bx[i] + AI_yy*m_By[i] + AI_yz*m_Bz[i]);
      grad[2] = inv_det*(AI_zx*m_Bx[i] + AI_zy*m_By[i] + AI_zz*m_Bz[i]);
    }
  }
};

} // namespace EDL_NAMESPACE

// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

#include <random>

TEST_CASE("TLSqGrad3D with linear gradient (scalar function)")
{
  using namespace EDL_NAMESPACE;
  typedef float real;
  //
  int num_iter = 100;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist1(-3.0, 3.0);
  std::uniform_real_distribution<> dist2(-0.5, 0.5);
  //
  for (int iter = 0; iter < num_iter; ++iter) {
    real a  = dist1(gen);
    real b  = dist1(gen);
    real c  = dist1(gen);
    real d  = dist1(gen);
    real x0 = dist2(gen);
    real y0 = dist2(gen);
    real z0 = dist2(gen);
    real f0 = a*x0 + b*y0 + c*z0 + d;
    TLSqGrad3D<real,1> lsq(x0, y0, z0, f0);
    //
    real x[6] = { 1, -1,  0,  0,  0,  0};
    real y[6] = { 0,  0,  1, -1,  0,  0};
    real z[6] = { 0,  0,  0,  0,  1, -1};
    for (int i = 0; i < 6; ++i) {
      x[i] += dist2(gen);
      y[i] += dist2(gen);
      z[i] += dist2(gen);
    }
    //
    for (int i = 0; i < 6; ++i) {
      real f = a*x[i] + b*y[i] + c*z[i] + d;
      lsq.addPoint(x[i], y[i], z[i], f);
    }
    //
    real grad[3];
    lsq.computeGrad(grad);
    //
    REQUIRE(grad[0] == doctest::Approx(a));
    REQUIRE(grad[1] == doctest::Approx(b));
    REQUIRE(grad[2] == doctest::Approx(c));    
  }
}

TEST_CASE("TLSqGrad3D with linear gradient (vector function)")
{
  using namespace EDL_NAMESPACE;
  typedef float real;
  //
  int num_iter = 100;
  const int NUM_VARS_TLSQ = 10;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist1(-3.0, 3.0);
  std::uniform_real_distribution<> dist2(-0.5, 0.5);
  //
  for (int iter = 0; iter < num_iter; ++iter) {
    real a[NUM_VARS_TLSQ], b[NUM_VARS_TLSQ], c[NUM_VARS_TLSQ], d[NUM_VARS_TLSQ], f0[NUM_VARS_TLSQ];
    for (int i = 0; i < NUM_VARS_TLSQ; ++i) {
      a[i] = dist1(gen);
      b[i] = dist1(gen);
      c[i] = dist1(gen);
      d[i] = dist1(gen);
    }
    real x0 = dist2(gen);
    real y0 = dist2(gen);
    real z0 = dist2(gen);
    for (int i = 0; i < NUM_VARS_TLSQ; ++i) {
      f0[i] = a[i]*x0 + b[i]*y0 + c[i]*z0 + d[i];
    }
    TLSqGrad3D<real,NUM_VARS_TLSQ> lsq(x0, y0, z0, f0);
    //
    real x[6] = { 1, -1,  0,  0,  0,  0};
    real y[6] = { 0,  0,  1, -1,  0,  0};
    real z[6] = { 0,  0,  0,  0,  1, -1};
    for (int i = 0; i < 6; ++i) {
      x[i] += dist2(gen);
      y[i] += dist2(gen);
      z[i] += dist2(gen);
    }
    //
    for (int i = 0; i < 6; ++i) {
      real f[NUM_VARS_TLSQ];
      for (int j = 0; j < NUM_VARS_TLSQ; ++j) {
        f[j] = a[j]*x[i] + b[j]*y[i] + c[j]*z[i] + d[j];
      }
      lsq.addPoint(x[i], y[i], z[i], f);
    }
    //
    real grad[NUM_VARS_TLSQ][3];
    lsq.computeGrad(grad);
    //
    for (int i = 0; i < NUM_VARS_TLSQ; ++i) {
      REQUIRE(grad[i][0] == doctest::Approx(a[i]));
      REQUIRE(grad[i][1] == doctest::Approx(b[i]));
      REQUIRE(grad[i][2] == doctest::Approx(c[i]));    
    }
  }
}



#endif // TLSQGRAD_H