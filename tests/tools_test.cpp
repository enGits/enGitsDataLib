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
#include <cmath>
#include <random>

#include "edl/tools.h"
#include "edl/stringtools.h"
#include "edl/geometrytools.h"

int main()
{
  std::string txt = "This is, a , comma separated, string.";
  auto parts = edl::StringTools::split(txt, ",");
  for (auto part : parts) {
    std::cout << "\"" << part << "\"" << std::endl;
  }
  std::cout << edl::StringTools::toLower(txt) << std::endl;
  std::cout << edl::StringTools::toUpper(txt) << std::endl;
  //
  edl::real eps = 100*std::numeric_limits<edl::real>::epsilon();
  edl::vec3_t x(1,1,1);
  edl::vec3_t y(1+eps, 1+eps, 1+eps);
  auto result = edl::almostEqual(x, y);
  //
  return(0);
}
