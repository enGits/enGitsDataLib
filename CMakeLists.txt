# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                         +
# + ENGITS CONFIDENTIAL                                                     +
# + ===================                                                     +
# +                                                                         +
# + Copyright 2013-2014 enGits GmbH, Todtnau, Germany                       +
# + All Rights Reserved.                                                    +
# +                                                                         +
# + NOTICE: All information contained herein is, and remains the property   +
# +         of enGits GmbH and its suppliers, if any. The intellectual and  +
# +         technical concepts contained herein are proprietary to          +
# +         enGits GmbH and are protected by trade secret or copyright law. +
# +         Dissemination of this information or reproduction of this       +
# +         material is strictly forbidden unless prior written permission  +
# +         is obtained from enGits GmbH.                                   +
# +                                                                         +
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cmake_minimum_required(VERSION 2.8.11)

# limit configuration types (must be done before project() statement)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited config" FORCE)

project(enGitsDataLib)

set(A_DEBUG_BUILD OFF CACHE BOOL "Build a debug version?")
if (${A_DEBUG_BUILD})
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
  add_definitions(-DEDL_DEBUG)
else()
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

set(A_SINGLE_PRECISION CACHE BOOL "Use single precision floating point data?")
if (${A_SINGLE_PRECISION})
  add_definitions(-DREAL_FLOAT)
endif()

set(A_ENABLE_VAR3 ON CACHE BOOL "Enable var3_t")
if (${A_ENABLE_VAR3})
  add_definitions(-DWITH_VAR3)
endif()

option(USE_QT5 "Compile using the Qt5 library." ON)

option(CREATE_STATIC_LIBRARY "Create a static version of engitsdatalib." OFF)

set(CMAKE_CXX_STANDARD 11)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX $ENV{HOME}/local CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(A_VERBOSE_BUILD OFF CACHE BOOL "Verbose build?")
if (${A_VERBOSE_BUILD})
  set(CMAKE_VERBOSE_MAKEFILE on)
else()
  set(CMAKE_VERBOSE_MAKEFILE off)
endif()

set(ENGITSDATALIB_MAJOR_VERSION 0)
set(ENGITSDATALIB_MINOR_VERSION 1)
set(ENGITSDATALIB_PATCH_VERSION 0)
set(ENGITSDATALIB_VERSION ${ENGITSDATALIB_MAJOR_VERSION}.${ENGITSDATALIB_MINOR_VERSION}.${ENGITSDATALIB_PATCH_VERSION})

# If windows, use cmath
# If Visual Studio, export all symbols for dlls and compile in parallel
if(WIN32)
  add_compile_definitions(_USE_MATH_DEFINES)
  if(MSVC)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    set(BUILD_SHARED_LIBS ON)
    add_compile_options(/MP)
  endif(MSVC)
endif(WIN32)

# set include directory, add src directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})  # This, according to the new CMake, shouldn't be needed...

# Qt4/Qt5
if(USE_QT5)
  add_definitions(-DWITH_QT5)
  set(CMAKE_AUTORCC ON)
  find_package(Qt5 COMPONENTS Core REQUIRED)
else()
  find_package(Qt4 COMPONENTS QtCore REQUIRED)
  include(${QT_USE_FILE})
  include_directories(${QT_QTCORE_INCLUDE_DIR})
endif(USE_QT5)

set(engitsdatalib_HEADERS
edl/bmap.h
edl/containertricks.h
edl/csvreader.h
edl/edlerror.h
edl/edl.h
edl/geometrytools.h
edl/interpolation.h
edl/interpolationcurve.h
edl/linsolve.h
edl/list.h
edl/mathvector.h
edl/mathvector_methods.h
edl/mathvector_operators.h
edl/mathvector_structs.h
edl/pmap.h
edl/pointcloudsurface.h
edl/profile.h
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
edl/tgraph.h
edl/tools.h
)

set(engitsdatalib_SOURCES
${engitsdatalib_HEADERS}
edl/csvreader.cpp
edl/edlerror.cpp
edl/geometrytools.cpp
edl/interpolationcurve.cpp
edl/list.cpp
edl/profile.cpp
edl/sparsetwodimarray.cpp
edl/timewindow.cpp
edl/unit.cpp
edl/varcontainer.cpp
)

if(CREATE_STATIC_LIBRARY)
  add_library(engitsdatalib STATIC ${engitsdatalib_SOURCES})
else(CREATE_STATIC_LIBRARY)
  add_library(engitsdatalib SHARED ${engitsdatalib_SOURCES})
endif(CREATE_STATIC_LIBRARY)

# If Visual Studio, generate export header to allow sharing of dll data..
if(WIN32 AND MSVC)
  include(GenerateExportHeader)
  generate_export_header(engitsdatalib)
endif(WIN32 AND MSVC)

set_target_properties(engitsdatalib PROPERTIES MACOSX_RPATH 1)
set_target_properties(engitsdatalib PROPERTIES VERSION ${ENGITSDATALIB_VERSION})
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set_target_properties(engitsdatalib PROPERTIES COMPILE_FLAGS "-Wno-deprecated")
endif()

if(USE_QT5)
  target_link_libraries(engitsdatalib Qt5::Core)
else()
  target_link_libraries(engitsdatalib ${QT_LIBRARIES})
  target_link_libraries(engitsdatalib Qt4::QtCore)
endif(USE_QT5)

install(
  TARGETS engitsdatalib
  LIBRARY DESTINATION lib${LIB_SUFFIX}
  ARCHIVE DESTINATION lib${LIB_SUFFIX}
  RUNTIME DESTINATION bin
)
install(FILES ${engitsdatalib_HEADERS} DESTINATION include/edl)
if(WIN32 AND MSVC)
  install(FILES ${PROJECT_BINARY_DIR}/engitsdatalib_export.h DESTINATION include/edl)
endif(WIN32 AND MSVC)

add_subdirectory(tests)
