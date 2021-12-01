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

#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>

#include "edl/lsqinterpolation.h"

using namespace edl;

typedef float test_t;
typedef MathVector<StaticVector<test_t,3>> vec_t;

test_t testfunc(vec_t x)
{
  test_t a0 = 1.0;
  test_t a1 = 1.5;
  test_t a2 = 2.0;
  test_t a3 = 2.5;
  return a0 + a1*sqr(x[0]) + a2*x[1] + a3*x[2];
}

int main()
{
  test_t h = 0.025;
  LsqInterpolation<test_t,3> lsq;
  vec_t X0(0*h, 0, 0);
  std::vector<vec_t> X;
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      for (int k = -1; k <= 1; ++k) {
        vec_t x(i*h,j*h,k*h);
        //std::cout << x << std::endl;
        X.push_back(x);
      }
    }
  }
  //
  std::vector<real> values;
  for (auto x : X) {
    test_t d = std::max(h/10000, (x - X0).abs());
    test_t w = 1.0/d;
    lsq.addNode(x, w);
    values.push_back(0);
  }
  auto   weights      = lsq.getMetrics(X0);
  test_t total_weight = 0;
  for (int i = 0; i < lsq.size(); ++i) {
    lsq.setValue(i, testfunc(X[i]));
    values[i] = testfunc(X[i]);
    total_weight += weights[i];
    std::cout << weights[i] << std::endl;
  }


  std::ofstream f("lsq.csv");
  lsq.write(f);
  std::cout << "\ncoefficients: " << lsq.coeffs() << std::endl;
  std::cout << weights.size() << " cells" << std::endl;
  std::cout << "total weight: " << total_weight << std::endl;
  std::cout << lsq.interpolate(X0) << std::endl;
  std::cout << testfunc(X0) << std::endl;

  real C_value = 0;
  for (int i = 0; i < lsq.size(); ++i) {
    C_value += weights[i]*values[i];
  }
  std::cout << C_value << std::endl;

  //
  return(0);
}
