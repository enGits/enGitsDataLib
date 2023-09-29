// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2015-2023 enGits GmbH                                    +
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
#ifndef MULTIDIMVECTOR_H
#define MULTIDIMVECTOR_H

#include "edl/edl.h"
#include "edl/edlerror.h"

#include <cstddef>
#include <vector>
#include <bits/stdint-uintn.h>

namespace EDL_NAMESPACE
{

template <uint8_t DIM, class T>
class MultiDimVector
{

  std::vector<MultiDimVector<DIM-1, T>> m_Data;

public:

  void resize(const std::vector<size_t>& sizes, int dim = 0)
  {
    if (sizes.size() - dim != DIM) {
      throw EdlError("MultiDimVector::resize(): sizes.size() - dim != DIM");
    }
    m_Data.resize(sizes[dim]);
    for (size_t i = 0; i < sizes[dim]; ++i) {
      m_Data[i].resize(sizes, dim + 1);
    }
  }

  MultiDimVector<DIM-1, T>& operator[](const size_t i) { return m_Data[i]; }

  size_t size(const int dim) const 
  { 
    if (dim == 0) {
      return m_Data.size();
    } else {
      return m_Data[0].size(dim - 1);
    }
  }
};

template <class T>
class MultiDimVector<1, T> : public std::vector<T>
{
public:

  void resize(const std::vector<size_t>& sizes, int dim = 0)
  {
    if (sizes.size() - dim != 1) {
      throw EdlError("MultiDimVector::resize(): sizes.size() - dim != 1");
    }
    std::vector<T>::resize(sizes[dim]);
  }

  size_t size(const int dim) const 
  { 
    if (dim == 0) {
      return std::vector<T>::size();
    } else {
      throw EdlError("MultiDimVector::size(): dim != 0");
    }
  }
};


} // namespace

// ----------------------------------------------------------------------------
// TESTS
// ----------------------------------------------------------------------------

TEST_CASE("MultiDimVector")
{
  typedef float real;
  //
  edl::MultiDimVector<3, real> v1;
  v1.resize({20, 30, 40});
  for (int i = 0; i < v1.size(0); ++i) {
    for (int j = 0; j < v1.size(1); ++j) {
      for (int k = 0; k < v1.size(2); ++k) {
        v1[i][j][k] = 100.0*i + 10.0*j + 1.0*k;
        v1[i][j][k] *= v1[i][j][k];
      }
    }
  }
  for (int i = 0; i < v1.size(0); ++i) {
    for (int j = 0; j < v1.size(1); ++j) {
      for (int k = 0; k < v1.size(2); ++k) {
        real x = 100.0*i + 10.0*j + 1.0*k;
        CHECK(v1[i][j][k] == doctest::Approx(x*x));
      }
    }
  }
}

#endif // MULTIDIMVECTOR_H