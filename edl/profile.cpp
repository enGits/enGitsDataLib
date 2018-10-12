// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2015 enGits GmbH                                         +
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
