cmake_minimum_required (VERSION 3.1)
# Fix behavior of CMAKE_CXX_STANDARD when targeting macOS.
if (POLICY CMP0025)
  cmake_policy(SET CMP0025 NEW)
endif ()
project (Kdict)
# The version number.
set (Kdict_VERSION_MAJOR 0)
set (Kdict_VERSION_MINOR 1)
set(CMAKE_BUILD_TYPE Debug)

# Add inc directory
include_directories(inc)
#include_directories(libs/co2/include)
#include_directories(libs/libaco-1.2.4)

# Add src files
file(GLOB SOURCES "src/*.cc")
#message("before ${SOURCES}")
#string(APPEND SOURCES ";libs/libaco-1.2.4/aco.c")
#string(APPEND SOURCES ";libs/libaco-1.2.4/acosw.S")
#message("after ${SOURCES}")
#set_property(SOURCE "libs/libaco-1.2.4/acosw.S" PROPERTY LANGUAGE C)

# Python library
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# Boost
find_package(Boost COMPONENTS thread system context REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})
find_package(Threads)

# co2 library
#add_subdirectory(libs/co2)
#add_subdirectory(libs/libaco-1.2.4)


# set c++14
set (CMAKE_CXX_STANDARD 14)


add_executable(Kdict ${SOURCES})
target_link_libraries(Kdict
    ${PYTHON_LIBRARIES}
    ${Boost_LIBRARIES}
    )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")
        message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")

        add_custom_command(TARGET Kdict POST_BUILD
            COMMAND dsymutil Kdict
        )
    endif()
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()

# configure a header file to pass some of the CMake settings
# to the source code
#configure_file (
#  "${PROJECT_SOURCE_DIR}/KdictConfig.h.in"
#  "${PROJECT_BINARY_DIR}/KdictConfig.h"
#  )



