// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include <vector>
#include <iostream>
#include <fstream>

#include "edl/interpolation.h"
#include "edl/tlist.h"

typedef edl::Interpolation<edl::LinInterpFunc<edl::real>> interp_t;

int main()
{
  using namespace edl;
  using namespace std;
  //
  real x1 = -1;
  real x2 =  5;
  real dx = 0.1;
  //
  {
    TList<real> X(10,10), Y(10,10);
    X.newEntry() = 0;
    X.newEntry() = 1;
    X.newEntry() = 2;
    X.newEntry() = 3;
    X.newEntry() = 4;
    FORALL(i, X.) {
      Y.newEntry() = X[i]*X[i];
    }
    interp_t I(&X, &Y, true);
    {
      real x = x1;
      while (x <= x2) {
        cout << x << ", " << x*x << ", " << I.interpolate(x) << endl;
        x += dx;
      }
    }
  }
  //
  {
    vector<real> X = {0,1,2,3,4}, Y;
    for (auto x : X) Y.push_back(x*x);
    interp_t I;
    I.init(X, Y, true);
    {
      real x = x1;
      while (x <= x2) {
        cout << x << ", " << x*x << ", " << I.interpolate(x) << endl;
        x += dx;
      }
    }
  }
}

