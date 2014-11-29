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
cmake_minimum_required(VERSION 2.8)

# limit configuration types (must be done before project() statement)
set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited config" FORCE)

project(enGitsDataLib)

set(A_DEBUG_BUILD OFF CACHE BOOL "Build a debug version?")
if (${A_DEBUG_BUILD})
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
else()
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

set(A_SINGLE_PRECISION CACHE BOOL "Use single precision floating point data?")
if (${A_SINGLE_PRECISION})
  add_definitions(-DREAL_FLOAT)
endif()

if (NOT INITIAL_CONFIG)
  set(INITIAL_CONFIG "..." CACHE INTERNAL "hidden flag...")
  set(CMAKE_INSTALL_PREFIX $ENV{HOME}/local CACHE STRING "Install path prefix, prepended onto install directories." FORCE)
endif()

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

# set include directory, add src directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

SET(engitsdatalib_HEADERS
edl.h
edlerror.h
list.h
tlist.h
)

SET(engitsdatalib_SOURCES
${engitsdatalib_HEADERS}
edlerror.cpp
list.cpp
)

ADD_LIBRARY(engitsdatalib SHARED ${engitsdatalib_SOURCES})

SET_TARGET_PROPERTIES(engitsdatalib PROPERTIES MACOSX_RPATH 1)
SET_TARGET_PROPERTIES(engitsdatalib PROPERTIES VERSION ${ENGITSDATALIB_VERSION})
SET_TARGET_PROPERTIES(engitsdatalib PROPERTIES COMPILE_FLAGS "-Wno-deprecated")

install(TARGETS engitsdatalib LIBRARY DESTINATION lib)
install(FILES ${engitsdatalib_HEADERS} DESTINATION include/edl)
