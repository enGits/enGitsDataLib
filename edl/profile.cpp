// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
// +                                                                    +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "profile.h"

namespace EDL_NAMESPACE
{
#if __cplusplus >= 201103L
  Profile::clock_t::time_point edl::Profile::m_Time[10];
  Profile::clock_t::duration edl::Profile::m_Duration[10] = {
    Profile::clock_t::duration(),
    Profile::clock_t::duration(),
    Profile::clock_t::duration(),
    Profile::clock_t::duration(),
    Profile::clock_t::duration(),
    Profile::clock_t::duration(),
    Profile::clock_t::duration(),
    Profile::clock_t::duration(),
    Profile::clock_t::duration(),
    Profile::clock_t::duration()
  };
#endif
}
