cmake_minimum_required (VERSION 3.1)
include(ExternalProject)

# Fix behavior of CMAKE_CXX_STANDARD when targeting macOS.
if (POLICY CMP0025)
    cmake_policy(SET CMP0025 NEW)
endif ()

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release")
endif()

project(kcollections)
# The version number.
set (Kdict_VERSION_MAJOR 0)
set (Kdict_VERSION_MINOR 1)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(SOURCE_DIR "kcollections/src")

message("looking in ${CMAKE_MODULE_PATH}")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

message("Using Python version " ${PYTHON_VERSION})

# Add inc directories
include_directories(inc)
include_directories(libs/uint256_t)

# Add src files
string(APPEND SOURCES "${SOURCE_DIR}/MurmurHash3.cc;")
string(APPEND SOURCES "${SOURCE_DIR}/UContainer.cc;")
string(APPEND SOURCES "${SOURCE_DIR}/Vertex.cc;")
string(APPEND SOURCES "${SOURCE_DIR}/Kcontainer.cc;")
string(APPEND SOURCES "libs/uint256_t/uint128_t.cpp;")
string(APPEND SOURCES "libs/uint256_t/uint256_t.cpp;")

string(APPEND PYTHON_CXX_SOURCES "${SOURCE_DIR}/CContainer.cc;")
string(APPEND DICT_SOURCES "${SOURCE_DIR}/Kdict.cc;")
string(APPEND SET_SOURCES "${SOURCE_DIR}/Kset.cc;")
string(APPEND COUNTER_SOURCES "${SOURCE_DIR}/Kcounter.cc;")

export(PACKAGE kcollections)

add_library(kcollections ${SOURCES})

include_directories(${PYTHON_INCLUDE_DIR})

# set c++14
set (CMAKE_CXX_STANDARD 14)

option(PYTHON "compile the python bindings")
if (${PYTHON})
    # pybind11 library
    add_subdirectory(libs/pybind11-2.2.4)
    set(PYBIND11_PYTHON_VERSION ${PYTHON_VERSION})
    include_directories(libs/pybind11-2.2.4/include)

    pybind11_add_module( _Kset ${SOURCES} ${SET_SOURCES} ${PYTHON_CXX_SOURCES})
    set_target_properties( _Kset PROPERTIES COMPILE_DEFINITIONS "KSET" )

    pybind11_add_module( _Kdict ${SOURCES} ${DICT_SOURCES} ${PYTHON_CXX_SOURCES})
    set_target_properties( _Kdict PROPERTIES COMPILE_DEFINITIONS "KDICT" )
    pybind11_add_module( _Kcounter ${SOURCES} ${COUNTER_SOURCES} ${PYTHON_CXX_SOURCES})
    set_target_properties( _Kcounter PROPERTIES COMPILE_DEFINITIONS "KCOUNTER" )
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")
        message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
    endif()
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  string( APPEND CMAKE_CXX_FLAGS " -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free" )
  message( "Using CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

# Jemalloc
set(PROJECT_JEMALLOC_PATH ${CMAKE_BINARY_DIR}/jemalloc-5.1.0)
ExternalProject_Add(project_jemalloc
  URL https://github.com/jemalloc/jemalloc/releases/download/5.1.0/jemalloc-5.1.0.tar.bz2
  URL_HASH MD5=1f47a5aff2d323c317dfa4cf23be1ce4
  PREFIX ${PROJECT_JEMALLOC_PATH}
  CONFIGURE_COMMAND ${PROJECT_JEMALLOC_PATH}/src/project_jemalloc/configure --prefix=${PROJECT_JEMALLOC_PATH}
  BUILD_IN_SOURCE TRUE
  BUILD_COMMAND ${MAKE})
ExternalProject_Get_Property(project_jemalloc install_dir)
add_library(jemalloc STATIC IMPORTED)
add_library(jemalloc_shared SHARED IMPORTED)
set_property(TARGET jemalloc PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libjemalloc${CMAKE_STATIC_LIBRARY_SUFFIX})
set_property(TARGET jemalloc PROPERTY INTERFACE_LINK_LIBRARIES pthread dl)
set_property(TARGET jemalloc_shared PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libjemalloc${CMAKE_SHARED_LIBRARY_SUFFIX})
include_directories(test_set_driver "${install_dir}/include")
if (${PYTHON})
    include_directories(_Kdict "${install_dir}/include")
    include_directories(_Kset "${install_dir}/include")

    add_dependencies(_Kdict project_jemalloc)
    add_dependencies(_Kset project_jemalloc)
    add_dependencies(_Kcounter project_jemalloc)
endif()
