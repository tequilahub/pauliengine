cmake_minimum_required(VERSION 3.28...4.2.1)

if(SKBUILD)
  # this value doesn't really matter, as it will be set by scikit-build-core
  set(PROJECT_VERSION "0.1.0")
else()
  include(FetchContent)
  FetchContent_Declare(
    CMakeExtraUtils
    QUIET
    URL
      https://github.com/LecrisUT/CMakeExtraUtils/archive/refs/tags/v0.4.1.tar.gz
      SOURCE_DIR
      ${CMAKE_CURRENT_BINARY_DIR}/external/upstream/_srcs/CMakeExtraUtils
  )
  FetchContent_MakeAvailable(CMakeExtraUtils)

  include(DynamicVersion)
  dynamic_version(PROJECT_PREFIX pauliengine_)
endif()

project(pauliengine LANGUAGES CXX VERSION ${PROJECT_VERSION})

# if CMAKE_BUILD_TYPE undefined, we set it to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()


# Options handling utilities Macro for printing an option in a consistent manner
# Written by Lori A. Burns (@loriab) and Ryan M. Richard (@ryanmrichard) Syntax:
# print_option(<option to print> <was specified>)
macro(print_option variable default)
  if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
    message(STATUS "Setting (unspecified) option ${variable}: ${default}")
  else()
    message(STATUS "Setting option ${variable}: ${${variable}}")
  endif()
endmacro()

# Wraps an option with default ON/OFF. Adds nice messaging to option() Written
# by Lori A. Burns (@loriab) and Ryan M. Richard (@ryanmrichard) Syntax:
# option_with_print(<option name> <description> <default value>)
macro(option_with_print variable msge default)
  print_option(${variable} ${default})
  option(${variable} ${msge} ${default})
endmacro()

# Wraps an option with a default other than ON/OFF and prints it Written by Lori
# A. Burns (@loriab) and Ryan M. Richard (@ryanmrichard) NOTE: Can't combine
# with above b/c CMake handles ON/OFF options specially NOTE2: CMake variables
# are always defined so need to further check for if they are the NULL string.
# This is also why we need the force Syntax: option_with_default(<option name>
# <description> <default value>)
macro(option_with_default variable msge default)
  print_option(${variable} "${default}")
  if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
    set(${variable} "${default}" CACHE STRING ${msge} FORCE)
  endif()
endmacro()

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ standard")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "C++ standard required")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "C++ extensions")

include(${PROJECT_SOURCE_DIR}/cmake/compiler_flags/CXXFlags.cmake)

message(STATUS "${PROJECT_NAME} version: ${PROJECT_VERSION_FULL}")

# report on compiler flags in use
message(STATUS "Configuring a ${CMAKE_BUILD_TYPE} build")
string(TOUPPER ${CMAKE_BUILD_TYPE} _cmake_build_type_upper)

message(STATUS "Compiler flags for ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "   From environment       : ${CMAKE_CXX_FLAGS}")
set(
  _cmake_build_type_specific_flags
  "${CMAKE_CXX_FLAGS_${_cmake_build_type_upper}}"
)
message(
  STATUS
  "   Build-type-specific    : ${_cmake_build_type_specific_flags}"
)
message(STATUS "   Vectorization flag     : ${ARCH_FLAG}")
message(
  STATUS
  "   Project defaults       : ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION} ${pauliengine_CXX_FLAGS}"
)
message(STATUS "   User-appended          : ${EXTRA_CXXFLAGS}")

include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})

configure_file(
  ${PROJECT_SOURCE_DIR}/include/pauliengine/Info.h.in
  ${PROJECT_BINARY_DIR}/include/pauliengine/Info.h
)

# Conan Toolchain einbinden (nach: conan install .. --install-folder=build --build=missing)
include(${CMAKE_BINARY_DIR}/conan_toolchain.cmake OPTIONAL RESULT_VARIABLE _found)

add_subdirectory(src)

