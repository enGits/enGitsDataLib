# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                    +
# + This file is part of enGitsDataLib.                                +
# +                                                                    +
# + Copyright 2015 enGits GmbH                                         +
# +                                                                    +
# + enGitsDataLib is free software: you can redistribute it and/or     +
# + modify it under the terms of the GNU Lesser General Public License +
# + as published by the Free Software Foundation, either version 3 of  +
# + the License, or (at your option) any later version.                +
# +                                                                    +
# + enGitsDataLib is distributed in the hope that it will be useful,   +
# + but WITHOUT ANY WARRANTY; without even the implied warranty of     +
# + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      +
# + GNU Lesser General Public License for more details.                +
# +                                                                    +
# + You should have received a copy of the GNU Lesser General Public   +
# + License along with enGitsDataLib.                                  +
# + If not, see <http://www.gnu.org/licenses/>.                        +
# +                                                                    +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# add_executable(pointcloud_test pointcloud_test.cpp)
# target_link_libraries(pointcloud_test engitsdatalib)

# add_executable(tools_test tools_test.cpp)
# target_link_libraries(tools_test engitsdatalib)

# add_executable(csv_test csv_test.cpp)
# target_link_libraries(csv_test engitsdatalib)

# add_executable(lsq_test lsq_test.cpp)
# target_link_libraries(lsq_test engitsdatalib)

# add_executable(interp_test interp_test.cpp)
# target_link_libraries(interp_test engitsdatalib)

add_executable(geometrytools_test geometrytools_test.cpp)
target_link_libraries(geometrytools_test engitsdatalib)

add_executable(doctest_test doctest_test.cpp)
target_link_libraries(doctest_test engitsdatalib)
