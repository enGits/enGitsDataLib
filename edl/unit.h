// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef UNIT_H
#define UNIT_H

#include <vector>
#include <string>

#include "edl/edl.h"


namespace EDL_NAMESPACE
{


class Unit
{

protected: // attributes

  real             m_Factor;
  std::vector<int> m_Units;


public: // methods

  Unit(std::string name);

  bool unitMatch(const Unit& unit) const;
  real convert(real value, const Unit& unit) const;

};

}

#endif
