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
  return a0 + a1*x[0] + a2*x[1] + a3*x[2];
}

int main()
{
  LsqInterpolation<test_t,3> lsq;
  std::vector<vec_t> X;
  X.push_back(vec_t(1,0,0));
  X.push_back(vec_t(1,1,0));
  X.push_back(vec_t(1,1,0.1));
  X.push_back(vec_t(0,0,1));
  //
  for (auto x : X) {
    lsq.addNode(x);
  }
  for (int i = 0; i < lsq.size(); ++i) {
    lsq.setValue(i, testfunc(X[i]));
  }

  std::ofstream f("lsq.csv");
  lsq.write(f);
  std::cout << lsq.coeffs() << std::endl;
  //
  return(0);
}
