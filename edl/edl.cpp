// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "edl/edl.h"
#define DOCTEST_CONFIG_IMPLEMENT
#include "edl/doctest.h"

namespace EDL_NAMESPACE
{

std::map<std::string, uint64_t> edl::smartBreakPoint::all_hits;

} // namespace