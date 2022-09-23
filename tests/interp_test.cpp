// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2022 enGits GmbH                                         +
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

