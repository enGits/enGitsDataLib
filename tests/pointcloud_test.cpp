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

#include "edl/pointcloudsurface.h"
#include "edl/geometrytools.h"

using namespace std;
using namespace edl;

int main()
{
  random_device rd;
  mt19937 gen(rd());
  //
  vector<vec3_t> X, N, IX;
  //
  int    num_long = 36;
  int    num_lat  = 18;
  real   radius   = 1.0;
  vec3_t n0(1,0,0);
  real   n_scale = 0.1;
  //
  uniform_real_distribution<edl::real> random(0, 0.001*radius);
  //
  n0.normalise();
  for (int i = 0; i < num_long; ++i) {
    edl::real alpha = 2*i*M_PI/num_long;
    for (int j = 0; j <= num_lat; ++j) {
      edl::real beta = 0.99*(j*M_PI/num_lat - 0.5*M_PI);
      vec3_t x(0, 0, radius);
      x = rotate(x, vec3_t(1,0,0), beta);
      x = rotate(x, vec3_t(0,0,1), alpha);
      real xn = x*n0;
      x -= xn*n0;
      xn *= n_scale;
      x += xn*n0;
      x[0] += 2;
      X.push_back(x);
      N.push_back(x.normalised());
    }
  }  
  edl::PointCloudSurface<edl::real> PC;
  PC.setPoints(X);
  PC.planeFit();
  PC.setNormals(PC.planeNormal());
  cout << "n_plane = " << PC.planeNormal() << endl;

  for (vec3_t x : X) {
    IX.push_back(vec3_t(x[0] + random(gen), x[1] + random(gen), x[2] + random(gen)));
  }
  
  PC.setInterpolationPoints(IX, PC.planeNormal());
  PC.computeWeights<PointCloudSurface<edl::real>::DistanceWeightedLeastSquares>();

  vector<edl::real> f_src, f_dst(IX.size());
  for (vec3_t x : X) {
    f_src.push_back(x[2]);
  }
  PC.interpolate(f_src, f_dst);
  real max_err = 0;
  for (int i = 0; i < X.size(); ++i) {
    max_err = max(max_err, fabs(f_src[i] - f_dst[i]));
  }
  cout << "max_err = " << max_err << endl;

  return(0); 
}
