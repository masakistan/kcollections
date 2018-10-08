cmake_minimum_required (VERSION 3.1)
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
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

project (Kdict)
# The version number.
set (Kdict_VERSION_MAJOR 0)
set (Kdict_VERSION_MINOR 1)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

message( "looking in ${CMAKE_MODULE_PATH}" )

# Add inc directory
include_directories(inc)
include_directories(libs/pybind11-2.2.4/include)

# Add src files
#file(GLOB SOURCES "src/*.cc")
#set(SOURCES "src/Bkmer.cc;")
#string(APPEND SOURCES "src/BloomFilter.cc;")
string(APPEND SOURCES "src/CContainer.cc;")
#string(APPEND SOURCES "src/Container.cc;")
#string(APPEND SOURCES "src/Kdict.cc;")
string(APPEND SOURCES "src/MurmurHash3.cc;")
#string(APPEND SOURCES "src/SufClustData.cc;")
string(APPEND SOURCES "src/UContainer.cc;")
string(APPEND SOURCES "src/Vertex.cc;")

# Python library
find_package(PythonLibs 2.7 REQUIRED)
find_package(PythonInterp 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIR})

# Boost
#find_package(Boost COMPONENTS thread system context chrono date_time atomic REQUIRED)
#include_directories(${Boost_INCLUDE_DIR})
#message("boost ${Boost_INCLUDE_DIR}")

# pybind11 library
add_subdirectory(libs/pybind11-2.2.4)


# set c++14
set (CMAKE_CXX_STANDARD 14)


#add_executable(check "src/TestDriver.cc;" ${SOURCES})
#target_link_libraries(check
#    ${PYTHON_LIBRARIES}
#    #${Boost_LIBRARIES}
#    )

#add_executable(test_driver "src/TestExternalFile.cc;" ${SOURCES})
#target_link_libraries(test_driver
#    ${PYTHON_LIBRARIES}
#    #${Boost_LIBRARIES}
#    )

add_executable( test_driver "src/TestLosePoint.cc;" ${SOURCES})



set(PYBIND11_PYTHON_VERSION 2.7)
#set(bindings_python_version 2.7)
#pybind11_add_module(kc ${SOURCES} ${Boost_LIBRARIES})

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")
        message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")

        add_custom_command(TARGET test_driver POST_BUILD
            COMMAND dsymutil test_driver
            )

        #add_custom_command(TARGET check POST_BUILD
        #    COMMAND dsymutil check
        #    )
    endif()
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()


