// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2015-2020 enGits GmbH                                    +
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

#ifndef TOOLS_H
#define TOOLS_H

#include "edl.h"
#include "edl/mathvector.h"

#include <cmath>
#include <cstddef>
#include <vector>

#include "doctest.h"

namespace EDL_NAMESPACE
{

template <typename T>
T interpolate1(T x, T x1, T x2, T y1, T y2)
{
  T w = (x - x1)/(x2 - x1);
  return (T(1) - w)*y1 + w*y2;
}

template <typename T>
T interpolate3(T x, T x1, T x2, T x3, T x4, T y1, T y2, T y3, T y4)
{
  T y_x0 = (x3-x2)*(y3 - y1)/(x3 - x1);
  T y_x1 = (x3-x2)*(y4 - y2)/(x4 - x2);
  //
  T a =                y_x0;
  T b =  3*(y3-y2) - 2*y_x0 - y_x1;
  T c = -2*(y3-y2) +   y_x0 + y_x1;
  //
  T w = (x - x2)/(x3 - x2);
  return a*w + b*w*w + c*w*w*w + y2;
}

template <class T, class C1, class C2>
T interpolate1D(T x, const C1& x_values, const C2& y_values) {
  // Check that x_values and y_values have the same size
  if (x_values.size() != y_values.size()) {
    throw std::invalid_argument("x_values and y_values must have the same size");
  }

  // Check that there are at least two points for interpolation
  if (x_values.size() < 2) {
    throw std::invalid_argument("At least two points are required for interpolation");
  }

  // Find the interval [x0, x1] that contains x
  auto it = std::lower_bound(x_values.begin(), x_values.end(), x);

  if (it == x_values.begin()) {
    // x is before the first element, use the first interval
    auto x0 = *it;
    auto x1 = *(it + 1);
    auto y0 = y_values[0];
    auto y1 = y_values[1];
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
  } else if (it == x_values.end()) {
    // x is beyond the last element, use the last interval
    auto x0 = *(it - 2);
    auto x1 = *(it - 1);
    auto y0 = y_values[y_values.size() - 2];
    auto y1 = y_values[y_values.size() - 1];
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
  } else {
    // x is between two elements, use the found interval
    auto idx = std::distance(x_values.begin(), it);
    auto x0 = x_values[idx - 1];
    auto x1 = x_values[idx];
    auto y0 = y_values[idx - 1];
    auto y1 = y_values[idx];
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
  }
}

template <typename C>
edl::StaticVector<typename C::value_type, 4> statistics(const C& container)
{
  typedef typename C::value_type scalar_t;
  //
  scalar_t sum(0);
  scalar_t sum_sq(0);
  scalar_t min_value =  std::numeric_limits<scalar_t>::max();
  scalar_t max_value = -std::numeric_limits<scalar_t>::max();
  for (auto i = container.begin(); i != container.end(); ++i) {
    auto value = *i;
    sum    += value;
    sum_sq += sqr(value);
    min_value = std::min(min_value, value);
    max_value = std::max(max_value, value);
  }
  //
  edl::StaticVector<typename C::value_type, 4> result;
  result[0] = sum/container.size();
  result[1] = std::sqrt(std::max(scalar_t(0), sum_sq/container.size() - sqr(result[0])));
  result[2] = min_value;
  result[3] = max_value;
  return result;
}

template <typename C> 
std::tuple<std::vector<typename C::value_type>, std::vector<typename C::value_type>> runningStatistics(const C& container, size_t ws)
{
  typedef typename C::value_type scalar_t;
  using namespace std;
  //
  size_t N = container.size();
  size_t n = N - ws + 1;
  if (n > N) {
    throw std::invalid_argument("runningStatistics: n must be less than the size of the container");
  }
  vector<scalar_t> mean(n, 0);
  vector<scalar_t> std_dev(n, 0);
  scalar_t sum    = 0;
  scalar_t sum_sq = 0;
  //
  for (size_t i = 0; i < N; ++i) {
    sum    += container[i];
    sum_sq += sqr(container[i]);
    if (i >= ws-1) {
      size_t j = i - ws + 1;
      mean[j]    = sum/ws;
      std_dev[j] = sqrt(max(scalar_t(0), sum_sq/ws - sqr(mean[j])));
      sum    -= container[j];
      sum_sq -= sqr(container[j]);
    }
  }
  //
  return make_tuple(mean, std_dev);
}

} // namespace



TEST_CASE("simple sine wave statistics")
{
  using namespace EDL_NAMESPACE;
  //
  std::vector<double> sine_wave;
  int N = 1000;
  for (int i = 0; i < N; ++i) {
    sine_wave.push_back(std::sin(i*2*M_PI/N));
  }
  auto stats = statistics(sine_wave);
  CHECK(stats[0] == doctest::Approx(0.0));
  CHECK(stats[1] == doctest::Approx(1.0/std::sqrt(2.0)));
  CHECK(stats[2] == doctest::Approx(-1.0));
  CHECK(stats[3] == doctest::Approx(1.0));
}

TEST_CASE("1D_interpolation")
{
  using namespace EDL_NAMESPACE;
  //
  std::vector<double> x_values = {1.0, 2.0, 3.0, 4.0};
  std::vector<double> y_values = {2.0, 3.0, 5.0, 4.0};

  CHECK(interpolate1D(0.0, x_values, y_values) == doctest::Approx(1.0));
  CHECK(interpolate1D(1.0, x_values, y_values) == doctest::Approx(2.0));
  CHECK(interpolate1D(2.5, x_values, y_values) == doctest::Approx(4.0));
  CHECK(interpolate1D(1.5, x_values, y_values) == doctest::Approx(2.5));
  CHECK(interpolate1D(5.0, x_values, y_values) == doctest::Approx(3.0));
}

TEST_CASE("simple running statistics")
{
  using namespace EDL_NAMESPACE;
  //
  std::vector<double> f;
  int N = 10;
  for (int i = 0; i < N; ++i) {
    f.push_back(1.0*i/(N-1));
  }
  auto [mean, std_dev] = runningStatistics(f, 3);
  for (size_t i = 0; i < mean.size(); ++i) {
    CHECK(mean[i]    == doctest::Approx(f[i+1]));
    CHECK(std_dev[i] == doctest::Approx(std_dev[0]));
  }
}


#endif // TOOLS_H
