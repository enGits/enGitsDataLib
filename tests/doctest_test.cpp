// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// +                                                                    +
// + Copyright 2023 enGits GmbH                                         +
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
#include "edl/doctest.h"

#define EDL_DOCTEST

#include "edl/tools.h"
#include "edl/fastremovelist.h"
#include "edl/tlsqgrad.h"
#include "edl/geometrytools.h"
#include "edl/pointcloudsurface.h"
#include "edl/lsqinterpolation.h"
#include "edl/weightedset.h"
#include "edl/octree.h"
#include "edl/multidimvector.h"
#include "edl/shortvector.h"
#include "edl/amr.h"
#include "edl/searchtree.h"

int main(int argc, const char** argv)
{
  doctest::Context context(argc, argv);
  int test_result = context.run();
  return test_result;
}