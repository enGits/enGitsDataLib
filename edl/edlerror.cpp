// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "edlerror.h"
#include <cstdlib>
#include <typeinfo>

namespace EDL_NAMESPACE
{

EdlError::EdlError(std::string message)
{
  m_Message = message;
  setName();
}

EdlError::~EdlError()
{
}

std::string EdlError::name()
{
  return m_Name;
}

void EdlError::setName()
{
  m_Name = typeid(*this).name();
}

} // namespace





