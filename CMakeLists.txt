cmake_minimum_required(VERSION 3.18)
project(kcollections LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(KCOLLECTIONS_OPENMP "Enable OpenMP for optional parallel benchmarks" OFF)
if(KCOLLECTIONS_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    set(KC_OPENMP_LIBS OpenMP::OpenMP_CXX)
  endif()
endif()

find_package(pybind11 CONFIG REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-omit-frame-pointer")
endif()

set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/kcollections/src")

pybind11_add_module(_kcollections MODULE
  "${SOURCE_DIR}/globals.cpp"
  "${SOURCE_DIR}/core/kset_core.cc"
  "${SOURCE_DIR}/core/kdict_core.cc"
  "${SOURCE_DIR}/core/kcounter_core.cc"
  "${SOURCE_DIR}/Kset.cc"
  "${SOURCE_DIR}/Kcounter.cc"
  "${SOURCE_DIR}/bindings_main.cc"
  "${SOURCE_DIR}/bindings_kset.cc"
  "${SOURCE_DIR}/bindings_kdict.cc"
  "${SOURCE_DIR}/bindings_kcounter.cc"
)

target_compile_definitions(_kcollections PRIVATE PYTHON)
target_link_libraries(_kcollections PRIVATE pybind11::pybind11 Threads::Threads ${KC_OPENMP_LIBS})
target_include_directories(_kcollections PRIVATE ${CMAKE_SOURCE_DIR}/inc)

install(TARGETS _kcollections LIBRARY DESTINATION kcollections)

option(KCOLLECTIONS_BUILD_CPP_TESTS "Build C++ unit test binaries" OFF)
if(KCOLLECTIONS_BUILD_CPP_TESTS)
  message(STATUS "C++ tests not configured for factored layout")
endif()
