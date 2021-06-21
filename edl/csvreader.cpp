// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2016 enGits GmbH                                         +
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

#include "csvreader.h"
#include "tools.h"
#include "edlerror.h"

#include <fstream>

namespace EDL_NAMESPACE
{

CsvReader::CsvReader(std::string file_name)
{
  std::ifstream f(file_name);
  std::string line;
  std::getline(f, line);
  auto header_parts = edl::split(line, ",");
  //
  while (!f.eof() && line != "") {
    std::getline(f, line);
    auto parts = edl::split(line, ",");
    for (int i = 0; i < parts.size(); ++i) {
      m_Data[header_parts[i]].push_front(parts[i]);
    }
  }
  //
  int N = -1;
  for(auto i : m_Data) {
    int n = 0;
    for (auto j : i.second) {
      ++n;
    }
    if (N < 0) {
      N = n;
    } else {
      if (N != n) {
        throw EdlError("Error while reading CSV file.");
      }
    }
  }
  if (N < 0) {
    throw EdlError("Error while reading CSV file.");
  }
}

}

