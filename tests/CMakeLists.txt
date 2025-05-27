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

add_executable(doctest_test doctest_test.cpp)
target_link_libraries(doctest_test engitsdatalib)
