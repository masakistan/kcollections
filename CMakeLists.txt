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

message( "looking in ${CMAKE_MODULE_PATH}" )
string( APPEND CMAKE_MODULE_PATH "/usr/local/Cellar/python/2.7.14/Frameworks/Python.framework/Versions/2.7/lib/" )

# Add inc directory
include_directories(inc)
include_directories(libs/pybind11-2.2.4/include)

# Add src files
file(GLOB SOURCES "src/*.cc")

# Python library
find_package(PythonLibs 2.7 REQUIRED)
find_package(PythonInterp 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# Boost
find_package(Boost COMPONENTS thread system context REQUIRED)
include_directories(${BOOST_INCLUDE_DIRS})

# pybind11 library
add_subdirectory(libs/pybind11-2.2.4)


# set c++14
set (CMAKE_CXX_STANDARD 14)


add_executable(Kdict ${SOURCES})
target_link_libraries(Kdict
    ${PYTHON_LIBRARIES}
    ${Boost_LIBRARIES}
    )

set(PYBIND11_PYTHON_VERSION 2.7)
#set(bindings_python_version 2.7)
pybind11_add_module(kc ${SOURCES} ${Boost_LIBRARIES})

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


