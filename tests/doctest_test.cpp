// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                    +
// + This file is part of enGitsDataLib.                                +
// + Copyright 2015-2025 enGits GmbH                                    +
// +                                                                    +
// + enGitsDataLib is released under the MIT License.                   +
// + See LICENSE file for details.                                      +
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
#include "edl/optimise.h"

int main(int argc, const char** argv)
{
  doctest::Context context(argc, argv);
  int test_result = context.run();
  return test_result;
}