# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                    +
# + This file is part of enGitsDataLib.                                +
# + Copyright 2015-2026 enGits GmbH                                    +
# +                                                                    +
# + enGitsDataLib is released under the MIT License.                   +
# + See LICENSE file for details.                                      +
# +                                                                    +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cmake_minimum_required(VERSION 3.22)

project(enGitsDataLib VERSION 1.0.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)
include(TestBigEndian)

# export build setup for development environment
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(A_DEBUG_BUILD "Build a debug version (deprecated; use CMAKE_BUILD_TYPE)." OFF)
if(A_DEBUG_BUILD AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
elseif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build.")
endif()
if(NOT CMAKE_CONFIGURATION_TYPES)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

option(A_SINGLE_PRECISION "Use single precision floating point data?" OFF)
option(A_ENABLE_VAR3 "Enable var3_t" ON)
option(EDL_ENABLE_EXPENSIVE_TESTS "Enable resource-intensive tests." OFF)

option(CREATE_STATIC_LIBRARY "Create a static version of engitsdatalib." OFF)

test_big_endian(PLATFORM_IS_BIG_ENDIAN)
if(PLATFORM_IS_BIG_ENDIAN)
  message(STATUS "System Endianness: big")
else()
  message(STATUS "System Endianness: little")
endif()

set(A_VERBOSE_BUILD OFF CACHE BOOL "Verbose build?")
if (${A_VERBOSE_BUILD})
  set(CMAKE_VERBOSE_MAKEFILE on)
else()
  set(CMAKE_VERBOSE_MAKEFILE off)
endif()

set(ENGITSDATALIB_MAJOR_VERSION ${PROJECT_VERSION_MAJOR})
set(ENGITSDATALIB_MINOR_VERSION ${PROJECT_VERSION_MINOR})
set(ENGITSDATALIB_PATCH_VERSION ${PROJECT_VERSION_PATCH})
set(ENGITSDATALIB_VERSION ${PROJECT_VERSION})

# If Visual Studio, export all symbols for dlls and compile in parallel
if(MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
  add_compile_options(/MP)
endif(MSVC)

set(engitsdatalib_HEADERS
edl/amr.h
edl/bmap.h
edl/containertricks.h
edl/csvreader.h
edl/doctest.h
edl/edlerror.h
edl/edl.h
edl/fastremovelist.h
edl/geometrytools.h
edl/interpolation.h
edl/interpolationcurve.h
edl/linsolve.h
edl/list.h
edl/lsqinterpolation.h
edl/lsqgrad.h
edl/mathvector.h
edl/mathvector_methods.h
edl/mathvector_operators.h
edl/mathvector_structs.h
edl/multidimvector.h
edl/octree.h
edl/optimise.h
edl/pmap.h
edl/pointcloudsurface.h
edl/profile.h
edl/quicktimer.h
edl/searchtree.h
edl/shortdeltavector.h
edl/shortptrvector.h
edl/shortvector.h
edl/smallsquarematrix.h
edl/sparsetwodimarray.h
edl/stringtools.h
edl/tkeylist.h
edl/tlist2d.h
edl/tlist.h
edl/timewindow.h
edl/tmappedvar.h
edl/tmdimindex.h
edl/tmdimlist.h
edl/tmdimvar.h
edl/tqueue.h
edl/tsdimvar.h
edl/tsparsetwodimarray.h
edl/unit.h
edl/varcontainer.h
edl/teelogger.h
edl/tgraph.h
edl/tlsqgrad.h
edl/tools.h
edl/weightedset.h
)

set(engitsdatalib_SOURCES
${engitsdatalib_HEADERS}
edl/csvreader.cpp
edl/edl.cpp
edl/edlerror.cpp
edl/interpolationcurve.cpp
edl/list.cpp
edl/profile.cpp
edl/sparsetwodimarray.cpp
edl/teelogger.cpp
edl/timewindow.cpp
edl/unit.cpp
edl/varcontainer.cpp
)

if(CREATE_STATIC_LIBRARY)
  add_library(engitsdatalib STATIC ${engitsdatalib_SOURCES})
else(CREATE_STATIC_LIBRARY)
  add_library(engitsdatalib SHARED ${engitsdatalib_SOURCES})
endif(CREATE_STATIC_LIBRARY)

target_compile_features(engitsdatalib PUBLIC cxx_std_17)

target_compile_definitions(engitsdatalib
  PUBLIC
    $<$<PLATFORM_ID:Windows>:_USE_MATH_DEFINES>
    $<$<CONFIG:Debug>:EDL_DEBUG>
    $<$<BOOL:${A_SINGLE_PRECISION}>:REAL_FLOAT>
    $<$<BOOL:${A_ENABLE_VAR3}>:WITH_VAR3>
    $<$<BOOL:${PLATFORM_IS_BIG_ENDIAN}>:PLATFORM_BIG_ENDIAN>
    $<$<NOT:$<BOOL:${PLATFORM_IS_BIG_ENDIAN}>>:PLATFORM_LITTLE_ENDIAN>
)

target_include_directories(engitsdatalib
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# If Visual Studio, generate export header to allow sharing of dll data..
if(WIN32 AND MSVC)
  include(GenerateExportHeader)
  generate_export_header(engitsdatalib)
  if(CREATE_STATIC_LIBRARY)
    target_compile_definitions(engitsdatalib PUBLIC ENGITSDATALIB_STATIC_DEFINE)
  endif()
endif(WIN32 AND MSVC)

set_target_properties(engitsdatalib PROPERTIES
  DEBUG_POSTFIX d
  MACOSX_RPATH 1
  SOVERSION ${ENGITSDATALIB_MAJOR_VERSION}
  VERSION ${ENGITSDATALIB_VERSION}
)
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set_target_properties(engitsdatalib PROPERTIES COMPILE_FLAGS "-Wno-deprecated")
endif()

install(
  TARGETS engitsdatalib
  EXPORT enGitsDataLibTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES ${engitsdatalib_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/edl)
if(WIN32 AND MSVC)
  install(FILES ${PROJECT_BINARY_DIR}/engitsdatalib_export.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/edl)
endif(WIN32 AND MSVC)

install(
  EXPORT enGitsDataLibTargets
  NAMESPACE enGitsDataLib::
  FILE enGitsDataLibTargets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/enGitsDataLib
)

configure_package_config_file(
  cmake/enGitsDataLibConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/enGitsDataLibConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/enGitsDataLib
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/enGitsDataLibConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/enGitsDataLibConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/enGitsDataLibConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/enGitsDataLib
)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
