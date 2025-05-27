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

#include "edl/csvreader.h"

int main()
{
  edl::CsvReader csv("test.csv");
  auto foo = csv.getColumn<std::string>("Foo");
  auto B   = csv.getColumn<edl::real>("B");
  auto bar = csv.getColumn<int>("bar");
  for (auto s : foo) std::cout << s << std::endl;
  return(0);
}
