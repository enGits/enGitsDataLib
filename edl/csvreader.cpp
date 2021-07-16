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
#include "stringtools.h"
#include "edlerror.h"

#include <fstream>

namespace EDL_NAMESPACE
{

CsvReader::CsvReader(std::string file_name)
{
  std::ifstream f(file_name);
  std::string line;
  std::getline(f, line);
  auto header_parts = StringTools::quotedSplit(line, ",");
  //
  while (!f.eof() && line != "") {
    line = StringTools::readQuotedLine(f, '"');
    auto parts = StringTools::quotedSplit(line, ",");
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
        /*
        std::vector<std::string> column_reversed(n);
        std::vector<std::string> column(n);
        std::copy(i.second.begin(), i.second.end(), column_reversed.begin());
        std::reverse_copy(column_reversed.begin(), column_reversed.end(), column.begin());
        */
        throw EdlError("Error while reading CSV file.");
      }
    }
  }
  if (N < 0) {
    throw EdlError("Error while reading CSV file.");
  }
}

}

