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

#include "edl/unit.h"
#include "edl/containertricks.h"
#include "edl/edlerror.h"

namespace EDL_NAMESPACE
{

Unit::Unit(std::string name)
{
  // mass, length, time, temperature
  //
  //
  // length
  //
  if (name == "m") {
    vlinit(m_Units) = 0, 1;
    m_Factor = 1.0;
  }
  if (name == "km") {
    vlinit(m_Units) = 0, 1;
    m_Factor = 1000;
  }
  //
  // volumetric flow rate
  //
  if (name == "m3/s") {
    vlinit(m_Units) = 0, 3, -1;
    m_Factor = 1.0;
  }
  if (name == "m3/min") {
    vlinit(m_Units) = 0, 3, -1;
    m_Factor = 1.0/60.0;
  }
  if (name == "m3/h") {
    vlinit(m_Units) = 0, 3, -1;
    m_Factor = 1.0/3600.0;
  }
  if (m_Units.size() == 0) {
    throw EdlError("unknown unit " + name);
  }
  while (m_Units.size() < 4) {
    m_Units.push_back(0);
  }
}

bool Unit::unitMatch(const Unit &unit) const
{
  if (m_Units.size() != unit.m_Units.size()) {
    return false;
  }
  for (size_t i = 0; i < m_Units.size(); ++i) {
    if (m_Units[i] != unit.m_Units[i]) {
      return false;
    }
  }
  return true;
}

real Unit::convert(real value, const Unit &unit) const
{
  if (!unitMatch(unit)) {
    throw EdlError("incompatible units");
  }
  return value*m_Factor/unit.m_Factor;
}


}

