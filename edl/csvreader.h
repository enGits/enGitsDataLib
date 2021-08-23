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

#ifndef CSVREADER_H
#define CSVREADER_H

#include "edl/edl.h"
#include "edl/stringtools.h"

#include <vector>
#include <string>
#include <map>
#include <forward_list>


namespace EDL_NAMESPACE
{


class CsvReader
{

protected: // attributes

  std::map<std::string, std::forward_list<std::string>> m_Data;


public: // methods

  CsvReader(std::string file_name);

  template <typename T>
  std::vector<T> getColumn(std::string name)
  {
    std::vector<T> result;
    if (m_Data.find(name) != m_Data.end()) {
      int N = 0;
      for (auto i : m_Data[name]) {
        ++N;
      }
      result.resize(N);
      int idx = 0;
      for (auto i : m_Data[name]) {
        T value;
        StringTools::stringTo(i, value);
        result[idx] = value;
        ++idx;
      }
    }
    std::reverse(result.begin(), result.end());
    return result;
  }

  std::vector<std::string> getColumnNames();


};

}

#endif
