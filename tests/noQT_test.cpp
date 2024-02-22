// tests for removing QT library from EDL

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "edl/doctest.h"

#include "edl/geometrytools.h"

using namespace EDL_NAMESPACE;


TEST_CASE("planeFit") {
  std::vector<vec3_t> points{{0, 0, 0}, {0, 1, 0}, {1, -1, 0}, {1, 2, 0}};
  vec3_t x0, n;
  planeFit(points, x0, n);

  CHECK(x0[0] == 0.5);
  CHECK(x0[1] == 0.5);
  CHECK(x0[2] == 0);
  CHECK(n[0] == 0);
  CHECK(n[1] == 0);
  CHECK(n[2] == 1);
}

TEST_CASE("polyNormal") {
  std::vector<vec3_t> points{{0, 0, 0}, {0, 1, 0}, {1, -1, 0}, {1, 2, 0}};
  vec3_t n = polyNormal(points);

  CHECK(n[0] == 0);
  CHECK(n[1] == 0);
  CHECK(n[2] == 2);
}
