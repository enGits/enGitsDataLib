// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "edl/csvreader.h"
#include "edl/stringtools.h"
#include "edl/edlerror.h"

#include <fstream>

namespace EDL_NAMESPACE
{

CsvReader::CsvReader(std::string file_name)
{
  std::ifstream f(file_name);
  if (!f.is_open()) {
    std::string error_message = "File \"" + file_name + "\" not found.";
    throw EdlError(error_message);
  }
  std::string line;
  //std::getline(f, line);
  line = StringTools::readQuotedLine(f, '"');
  auto header_parts = StringTools::quotedSplit(line, ",");
  //
  std::vector<bool> active_column(header_parts.size(), false);
  while (!f.eof() && line != "") {
    line = StringTools::readQuotedLine(f, '"');
    if (!edl::StringTools::trim(line).empty()) {
      auto parts = StringTools::quotedSplit(line, ",");
      for (int i = 0; i < parts.size(); ++i) {
        auto header = header_parts[i];
        if (!header.empty()) {
          active_column[i] = true;
          m_Data[header].push_front(StringTools::trim(parts[i]));
        }
      }
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

std::vector<std::string> CsvReader::getColumnNames()
{
  std::vector<std::string> column_names;
  for(auto i : m_Data) {
    column_names.push_back(i.first);
  }
  return column_names;
}

}

