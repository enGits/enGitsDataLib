// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
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

