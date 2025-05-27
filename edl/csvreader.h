// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
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
