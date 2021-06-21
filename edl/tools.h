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
#include <string>
#include <vector>

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

inline std::string trim(std::string s)
{
  std::string result = "";
  for (auto c : s) {
    if (!isspace(c) || result.size()) {
      result += c;
    }
  }
  while (result.size()) {
    if (isspace(result.back())) {
      result.pop_back();
    } else {
      break;
    }
  }
  return result;
}

inline std::vector<std::string> split(std::string s, std::string delimiter, bool trim_result=true)
{
  size_t pos_start = 0;
  size_t pos_end;
  size_t delim_len = delimiter.length();
  //
  std::string              token;
  std::vector<std::string> result;
  //
  while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
    token = s.substr (pos_start, pos_end - pos_start);
    pos_start = pos_end + delim_len;
    if (trim_result) {
      result.push_back(trim(token));
    } else {
      result.push_back(token);
    }
  }
  //
  if (trim_result) {
    result.push_back(trim(s.substr(pos_start)));
  } else {
    result.push_back(s.substr(pos_start));
  }
  return result;
}


} // namespace

#endif // TOOLS_H
