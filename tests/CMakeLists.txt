# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                    +
# + This file is part of enGitsDataLib.                                +
# + Copyright 2015-2025 enGits GmbH                                    +
# +                                                                    +
# + enGitsDataLib is released under the MIT License.                   +
# + See LICENSE file for details.                                      +
# +                                                                    +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

add_executable(geometrytools_test geometrytools_test.cpp)
target_link_libraries(geometrytools_test engitsdatalib)
add_test(NAME geometrytools_test COMMAND geometrytools_test)

add_executable(doctest_test doctest_test.cpp)
target_link_libraries(doctest_test engitsdatalib)
configure_file(sphere_test.mesh sphere_test.mesh COPYONLY)
add_test(
  NAME doctest_test
  COMMAND doctest_test --test-case-exclude=AMRMesh_double_pyramid_mockup
)

if(EDL_ENABLE_EXPENSIVE_TESTS)
  add_test(
    NAME doctest_expensive_test
    COMMAND doctest_test --test-case=AMRMesh_double_pyramid_mockup
  )
endif()
