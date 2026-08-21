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
set(edl_tests geometrytools_test)

add_executable(doctest_test doctest_test.cpp)
target_link_libraries(doctest_test engitsdatalib)
configure_file(sphere_test.mesh sphere_test.mesh COPYONLY)
add_test(
  NAME doctest_test
  COMMAND doctest_test --test-case-exclude=AMRMesh_double_pyramid_mockup
)
list(APPEND edl_tests doctest_test)

if(EDL_ENABLE_EXPENSIVE_TESTS)
  add_test(
    NAME doctest_expensive_test
    COMMAND doctest_test --test-case=AMRMesh_double_pyramid_mockup
  )
  list(APPEND edl_tests doctest_expensive_test)
endif()

if(WIN32 AND NOT CREATE_STATIC_LIBRARY)
  if(CMAKE_VERSION VERSION_LESS 3.22)
    message(FATAL_ERROR "Windows shared-library tests require CMake 3.22 or newer")
  endif()
  set_tests_properties(${edl_tests} PROPERTIES
    ENVIRONMENT_MODIFICATION
      "PATH=path_list_prepend:$<TARGET_FILE_DIR:engitsdatalib>"
  )
endif()
