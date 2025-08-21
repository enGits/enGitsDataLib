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

#ifndef OPTIMISE_H
#define OPTIMISE_H

#include <cmath>
#include <cstdint>
#include <iostream>
#include "edl/edl.h"
#include "edl/edlerror.h"
#include "edl/mathvector.h"
#include "edl/smallsquarematrix.h"

namespace EDL_NAMESPACE
{

template <typename TFunc, uint8_t DIM, typename TValue=real>
class Optimiser
{

  protected: // data types

  typedef TValue                                value_t;
  typedef MathVector<StaticVector<value_t,DIM>> vector_t;
  typedef SmallSquareMatrix<value_t,DIM>        matrix_t;
  typedef TFunc                                 func_t;


  protected: // attributes

  vector_t m_Delta1;
  vector_t m_Delta2;
  bool     m_UseNewton = true;


protected: // methods

  void computeDeltas(const vector_t& x)
  {
    const value_t eps = std::numeric_limits<value_t>::epsilon();
    const value_t c1  = std::pow(eps, value_t(1.0/3.0));
    const value_t c2  = std::pow(eps, value_t(1.0/4.0));
    for (int i = 0; i < DIM; ++i) {
      const value_t scale = std::max(std::abs(x[i]), value_t(1.0));
      m_Delta1[i] = c1 * scale; // delta for gradient
      m_Delta2[i] = c2 * scale; // delta for Jacobian
    }
  }

  vector_t optimiseNewton(const vector_t& x0, value_t tol, int max_iter)
  {
    vector_t x     = x0;
    bool     done  = false;
    int      iter  = 1;
    value_t  relax = 1.0;
    //
    while (!done) {
      if (iter > max_iter/2) {
        relax *= 0.95; // reduce step size
        if (relax < 0.01) {
          relax = 0.01; // minimum step size
        }
      }
      if (iter > max_iter) {
        done = true; // stop after max_iter
      }
      matrix_t J;
      func_t   f;
      vector_t nabla_f;
      computeDeltas(x);
      //
      // compute gradient at the position x and check for convergence
      //
      for (int i = 0; i < DIM; ++i) {
        vector_t x_ip = x;
        x_ip[i] += m_Delta1[i];
        vector_t x_im = x;
        x_im[i] -= m_Delta1[i];
        //
        value_t f_ip = f(x_ip);
        value_t f_im = f(x_im);
        nabla_f[i] = (f_ip - f_im) / (2 * m_Delta1[i]);
      }
      //
      // compute Jacobian at the position x
      //
      for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) {
          vector_t x_jp = x;
          x_jp[j] += m_Delta2[j];
          vector_t x_jm = x;
          x_jm[j] -= m_Delta2[j];
          //
          vector_t x_jp_ip = x_jp;
          x_jp_ip[i] += m_Delta2[i];
          vector_t x_jp_im = x_jp;
          x_jp_im[i] -= m_Delta2[i];
          //
          vector_t x_jm_ip = x_jm;
          x_jm_ip[i] += m_Delta2[i];
          vector_t x_jm_im = x_jm;
          x_jm_im[i] -= m_Delta2[i];
          //
          value_t f_jp_ip = f(x_jp_ip);
          value_t f_jp_im = f(x_jp_im);
          value_t f_jm_ip = f(x_jm_ip);
          value_t f_jm_im = f(x_jm_im);
          //
          value_t df_dxi_jp = (f_jp_ip - f_jp_im) / (2 * m_Delta2[i]);
          value_t df_dxi_jm = (f_jm_ip - f_jm_im) / (2 * m_Delta2[i]);
          //
          J[i][j] = (f_jp_ip - f_jp_im - f_jm_ip + f_jm_im) / (4 * m_Delta2[i] * m_Delta2[j]);
        }
      }
      auto JI = J.inverse();
      auto Dx = JI * nabla_f;
      if (Dx.abs() < tol) {
        done = true;
      }
      //
      x -= relax * Dx;
      ++iter;
    }
    //
    return x;
  }

  vector_t optimiseGradientDescent(const vector_t& x0, value_t tol, int max_iter)
  {
    EDL_BUG;            
    vector_t x = x0;
    computeDeltas(x);
    for (int iter = 0; iter < max_iter; ++iter) {
      func_t   f;
      vector_t nabla_f;
      //
      // compute gradient at the position x
      //
      for (int i = 0; i < DIM; ++i) {
        vector_t x_ip = x;
        x_ip[i] += m_Delta1[i];
        vector_t x_im = x;
        x_im[i] -= m_Delta1[i];
        //
        value_t f_ip = f(x_ip);
        value_t f_im = f(x_im);
        nabla_f[i] = (f_ip - f_im) / (2 * m_Delta1[i]);
      }
      //
      if (nabla_f.abs() < tol) {
        break; // convergence reached
      }
      //
      // update position
      //
      x -= nabla_f;
    }
    //
    return x;
  }



public: // methods

  vector_t optimise(const vector_t& x0=vector_t(0), value_t tol=1e-6, int max_iter=1000)
  {
    if (m_UseNewton) {
      return optimiseNewton(x0, tol, max_iter);
    }
  }
};

} // namespace EDL_NAMESPACE

template<typename TValue, uint8_t DIM, int X0=0, int Y0=0, int Z0=0> struct TestFunc 
{
  typedef edl::MathVector<edl::StaticVector<TValue,DIM>> vec_t;
  TValue operator()(const vec_t& x) const
  {
    return (edl::sqr(x[0] - TValue(X0)) + edl::sqr(x[1] - TValue(Y0)));
  }
};

// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

#include <random>
#include <type_traits>

static constexpr int   NUM_TRIALS = 1000; // <— your single place to change trial count

template<typename Real>
struct Tolerances {
  static constexpr Real eps =
    std::is_same_v<Real,float> ? Real(1e-4f) : Real(1e-8);
};

// ——————————————————— Helpers (type-generic) ———————————————————
template<typename Real>
struct Rng {
  std::mt19937_64 gen;
  explicit Rng(uint64_t seed = 0xC0FFEEULL) : gen(seed) {}
  Real uni(Real a, Real b) { std::uniform_real_distribution<Real> d(a,b); return d(gen); }
  Real nrm(Real s = Real(1)) { std::normal_distribution<Real> d(Real(0), s); return d(gen); }
};

template<int DIM, typename Real>
static edl::SmallSquareMatrix<Real,DIM>
make_spd(Rng<Real>& rng,
         const edl::StaticVector<Real,DIM>& diag_min,
         const edl::StaticVector<Real,DIM>& diag_max,
         Real offscale = Real(0.25))
{
  using mat = edl::SmallSquareMatrix<Real,DIM>;
  mat Q = mat::zero();

  for (int i = 0; i < DIM; ++i) {
    const Real d = rng.uni(diag_min[i], diag_max[i]);
    Q[i][i] = d;
  }
  for (int k = 0; k < 2; ++k) {
    edl::StaticVector<Real,DIM> u;
    for (int i = 0; i < DIM; ++i) u[i] = offscale * rng.nrm(Real(1));
    for (int i = 0; i < DIM; ++i)
      for (int j = 0; j < DIM; ++j)
        Q[i][j] += u[i] * u[j];
  }
  return Q;
}

template<int DIM, typename Real>
static edl::MathVector<edl::StaticVector<Real,DIM>>
rand_start(Rng<Real>& rng,
           const edl::MathVector<edl::StaticVector<Real,DIM>>& c,
           Real span_small, Real span_large)
{
  using vec = edl::MathVector<edl::StaticVector<Real,DIM>>;
  vec x0;
  for (int i = 0; i < DIM; ++i) {
    const Real span = (i % 2 == 0) ? span_large : span_small;
    x0[i] = c[i] + rng.uni(-span, span);
  }
  return x0;
}

// —————————————————— Functors (option B style) ——————————————————
template<int DIM, typename Real>
struct QuadND {
  using vec = edl::MathVector<edl::StaticVector<Real,DIM>>;
  using mat = edl::SmallSquareMatrix<Real,DIM>;
  static vec c;   // center (minimizer)
  static mat Q;   // symmetric PD
  Real operator()(const vec& x) const {
    vec d  = x - c;
    vec qd = Q * d;
    Real v = Real(0);
    for (int i = 0; i < DIM; ++i) v += d[i] * qd[i];
    return v;
  }
};
template<int DIM, typename Real> typename QuadND<DIM,Real>::vec QuadND<DIM,Real>::c;
template<int DIM, typename Real> typename QuadND<DIM,Real>::mat QuadND<DIM,Real>::Q;

template<int DIM, typename Real>
struct ExpBowlND {
  using vec = edl::MathVector<edl::StaticVector<Real,DIM>>;
  static vec c;                                   // center
  static edl::StaticVector<Real,DIM> w;           // positive weights
  Real operator()(const vec& x) const {
    Real acc = Real(0);
    for (int i = 0; i < DIM; ++i) {
      Real d = x[i] - c[i];
      acc += w[i] * d * d;
    }
    return std::exp(acc) - Real(1);
  }
};
template<int DIM, typename Real> typename ExpBowlND<DIM,Real>::vec ExpBowlND<DIM,Real>::c;
template<int DIM, typename Real> edl::StaticVector<Real,DIM> ExpBowlND<DIM,Real>::w;

template<int DIM, typename Real>
struct NastyHole {
  using vec = edl::MathVector<edl::StaticVector<Real,DIM>>;
  static vec  c;
  static Real rate;
  static Real offset;
  Real operator()(const vec& x) const {
    auto dx = x - c;
    Real r  = dx.abs();
    return offset + rate*std::pow(r,1.2);
  }
};
template<int DIM, typename Real> typename NastyHole<DIM,Real>::vec NastyHole<DIM,Real>::c;
template<int DIM, typename Real> Real NastyHole<DIM,Real>::rate;
template<int DIM, typename Real> Real NastyHole<DIM,Real>::offset;

// —————————————————— Test bodies (templated on Real) ——————————————————
template<typename Real>
static void run_case_2d_diag()
{
  using vec2 = edl::MathVector<edl::StaticVector<Real,2>>;
  Rng<Real> rng(1234);
  edl::Optimiser<QuadND<2,Real>,2,Real> opt;

  for (int t = 0; t < NUM_TRIALS; ++t) {
    QuadND<2,Real>::c = vec2(rng.uni(Real(-10), Real(10)),
                             rng.uni(Real(-1e3), Real(1e3)));

    QuadND<2,Real>::Q = edl::SmallSquareMatrix<Real,2>::zero();
    QuadND<2,Real>::Q[0][0] = rng.uni(Real(0.5), Real(5.0));
    QuadND<2,Real>::Q[1][1] = rng.uni(Real(2.0), Real(50.0));

    vec2 x0 = rand_start<2,Real>(rng, QuadND<2,Real>::c, Real(0.1), Real(5));

    auto x  = opt.optimise(x0, Tolerances<Real>::eps);

    CHECK(x[0] == doctest::Approx(QuadND<2,Real>::c[0]).epsilon(Tolerances<Real>::eps));
    CHECK(x[1] == doctest::Approx(QuadND<2,Real>::c[1]).epsilon(Tolerances<Real>::eps));
  }
}

template<typename Real>
static void run_case_2d_cross()
{
  using vec2 = edl::MathVector<edl::StaticVector<Real,2>>;
  Rng<Real> rng(5678);
  edl::Optimiser<QuadND<2,Real>,2,Real> opt;

  for (int t = 0; t < NUM_TRIALS; ++t) {
    QuadND<2,Real>::c = vec2(rng.uni(Real(-100), Real(100)),
                             rng.uni(Real(-1e3), Real(1e3)));

    edl::StaticVector<Real,2> dmin, dmax;
    dmin[0]=Real(1.0); dmax[0]=Real(4.0);
    dmin[1]=Real(0.5); dmax[1]=Real(3.0);
    QuadND<2,Real>::Q = make_spd<2,Real>(rng, dmin, dmax, Real(0.5));

    vec2 x0 = rand_start<2,Real>(rng, QuadND<2,Real>::c, Real(0.2), Real(20));
    auto x  = opt.optimise(x0, Tolerances<Real>::eps);

    CHECK(x[0] == doctest::Approx(QuadND<2,Real>::c[0]).epsilon(Tolerances<Real>::eps));
    CHECK(x[1] == doctest::Approx(QuadND<2,Real>::c[1]).epsilon(Tolerances<Real>::eps));
  }
}

template<typename Real>
static void run_case_3d_diag()
{
  using vec3 = edl::MathVector<edl::StaticVector<Real,3>>;
  Rng<Real> rng(2468);
  edl::Optimiser<QuadND<3,Real>,3,Real> opt;

  for (int t = 0; t < NUM_TRIALS; ++t) {
    QuadND<3,Real>::c = vec3(rng.uni(Real(-1e-3), Real(1e-3)),
                             rng.uni(Real(-10), Real(10)),
                             rng.uni(Real(100), Real(1000)));

    QuadND<3,Real>::Q = edl::SmallSquareMatrix<Real,3>::zero();
    QuadND<3,Real>::Q[0][0] = rng.uni(Real(0.5), Real(2.0));
    QuadND<3,Real>::Q[1][1] = rng.uni(Real(5.0), Real(20.0));
    QuadND<3,Real>::Q[2][2] = rng.uni(Real(0.01), Real(0.1));

    vec3 x0 = rand_start<3,Real>(rng, QuadND<3,Real>::c, Real(0.05), Real(200));
    auto x  = opt.optimise(x0, Tolerances<Real>::eps);

    CHECK(x[0] == doctest::Approx(QuadND<3,Real>::c[0]).epsilon(Tolerances<Real>::eps));
    CHECK(x[1] == doctest::Approx(QuadND<3,Real>::c[1]).epsilon(Tolerances<Real>::eps));
    CHECK(x[2] == doctest::Approx(QuadND<3,Real>::c[2]).epsilon(Tolerances<Real>::eps));
  }
}

template<typename Real>
static void run_case_3d_coupled()
{
  using vec3 = edl::MathVector<edl::StaticVector<Real,3>>;
  Rng<Real> rng(1357);
  edl::Optimiser<QuadND<3,Real>,3,Real> opt;

  for (int t = 0; t < NUM_TRIALS; ++t) {
    QuadND<3,Real>::c = vec3(rng.uni(Real(-5), Real(5)),
                             rng.uni(Real(500), Real(1500)),
                             rng.uni(Real(-1), Real(1)));

    edl::StaticVector<Real,3> dmin, dmax;
    dmin[0]=Real(1.0); dmax[0]=Real(5.0);
    dmin[1]=Real(0.5); dmax[1]=Real(3.0);
    dmin[2]=Real(2.0); dmax[2]=Real(8.0);
    QuadND<3,Real>::Q = make_spd<3,Real>(rng, dmin, dmax, Real(0.4));

    vec3 x0 = rand_start<3,Real>(rng, QuadND<3,Real>::c, Real(0.5), Real(50));
    auto x  = opt.optimise(x0, Tolerances<Real>::eps);

    CHECK(x[0] == doctest::Approx(QuadND<3,Real>::c[0]).epsilon(Tolerances<Real>::eps));
    CHECK(x[1] == doctest::Approx(QuadND<3,Real>::c[1]).epsilon(Tolerances<Real>::eps));
    CHECK(x[2] == doctest::Approx(QuadND<3,Real>::c[2]).epsilon(Tolerances<Real>::eps));
  }
}

template<typename Real>
static void run_case_2d_exp()
{
  using vec2 = edl::MathVector<edl::StaticVector<Real,2>>;
  Rng<Real> rng(4242);
  edl::Optimiser<ExpBowlND<2,Real>,2,Real> opt;

  for (int t = 0; t < NUM_TRIALS; ++t) {
    ExpBowlND<2,Real>::c[0] = rng.uni(Real(-5), Real(5));
    ExpBowlND<2,Real>::c[1] = rng.uni(Real(-50), Real(50));
    ExpBowlND<2,Real>::w[0] = rng.uni(Real(0.01), Real(0.2));
    ExpBowlND<2,Real>::w[1] = rng.uni(Real(0.001), Real(0.05));

    vec2 x0 = rand_start<2,Real>(rng, ExpBowlND<2,Real>::c, Real(0.5), Real(10));
    auto x  = opt.optimise(x0, Tolerances<Real>::eps);

    CHECK(x[0] == doctest::Approx(ExpBowlND<2,Real>::c[0]).epsilon(5*Tolerances<Real>::eps));
    CHECK(x[1] == doctest::Approx(ExpBowlND<2,Real>::c[1]).epsilon(5*Tolerances<Real>::eps));
  }
}

// —————————————————— Doctest entry points (float & double) ——————————————————
TEST_CASE("Optimise__float__2D_quadratic_diagonal__random")   { run_case_2d_diag<float>(); }
TEST_CASE("Optimise__float__2D_quadratic_cross__random")      { run_case_2d_cross<float>(); }
TEST_CASE("Optimise__float__3D_quadratic_diagonal__random")   { run_case_3d_diag<float>(); }
TEST_CASE("Optimise__float__3D_quadratic_coupled__random")    { run_case_3d_coupled<float>(); }
TEST_CASE("Optimise__float__2D_exponential_bowl__random")     { run_case_2d_exp<float>(); }

TEST_CASE("Optimise__double__2D_quadratic_diagonal__random")  { run_case_2d_diag<double>(); }
TEST_CASE("Optimise__double__2D_quadratic_cross__random")     { run_case_2d_cross<double>(); }
TEST_CASE("Optimise__double__3D_quadratic_diagonal__random")  { run_case_3d_diag<double>(); }
TEST_CASE("Optimise__double__3D_quadratic_coupled__random")   { run_case_3d_coupled<double>(); }
TEST_CASE("Optimise__double__2D_exponential_bowl__random")    { run_case_2d_exp<double>(); }

#endif // OPTIMISE_H