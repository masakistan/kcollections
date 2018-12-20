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

message( "cmake moduel path ${CMAKE_MODULE_PATH}" )

project (Kdict)
# The version number.
set (Kdict_VERSION_MAJOR 0)
set (Kdict_VERSION_MINOR 1)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(SOURCE_DIR "kcollections/src")
#set(ROARING_BUILD_STATIC ON)
set(ENABLE_ROARING_TESTS OFF)

message( "looking in ${CMAKE_MODULE_PATH}" )

# set the default version of Python to use at 2.7
# can be changed via command line, e.g. cmake -DPython3=ON ..
option(Python2 "Python2" ON)
option(Python3 "Python3" OFF)

if(Python2)
  set(PYTHON_VERSION 2.7)
endif()
if(Python3)
  set(PYTHON_VERSION 3.4)
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
message("Using Python version " ${PYTHON_VERSION})

# Add inc directory
include_directories(${PROJECT_SOURCE_DIR}/inc)
include_directories(${PROJECT_SOURCE_DIR}/libs/pybind11-2.2.4/include)
include_directories(${PROJECT_SOURCE_DIR}/libs/CRoaring/cpp)
include_directories(${PROJECT_SOURCE_DIR}/libs/CRoaring/include)

include_directories(${PROJECT_SOURCE_DIR}/libs/uint256_t)

# Add src files
string(APPEND SOURCES "${SOURCE_DIR}/MurmurHash3.cc;")
#string(APPEND SOURCES "${SOURCE_DIR}/CContainer.cc;")
string(APPEND SOURCES "${SOURCE_DIR}/UContainer.cc;")
string(APPEND SOURCES "${SOURCE_DIR}/Vertex.cc;")
string(APPEND SOURCES "${SOURCE_DIR}/Kcollections.cc;")
string(APPEND SOURCES "${SOURCE_DIR}/Kcontainer.cc;")
string(APPEND SOURCES "libs/uint256_t/uint128_t.cpp;")
string(APPEND SOURCES "libs/uint256_t/uint256_t.cpp;")

string(APPEND DICT_SOURCES "${SOURCE_DIR}/Kdict.cc;")
string(APPEND SET_SOURCES "${SOURCE_DIR}/Kset.cc;")
string(APPEND COUNTER_SOURCES "${SOURCE_DIR}/Kcounter.cc;")
string(APPEND COLOR_SOURCES "${SOURCE_DIR}/Kcolor.cc;")
#add_executable(settest "${SOURCE_DIR}/TestSetLoad.cxx;" ${SET_SOURCES} ${SOURCES})

#add_executable(settestserial "${SOURCE_DIR}/TestSetLoadSerial.cxx;" ${SET_SOURCES} ${SOURCES})
#add_executable(settestserialfast "${SOURCE_DIR}/TestSetLoadSerialFast.cxx;" ${SET_SOURCES} ${SOURCES})
#target_link_libraries(settest ${PYTHON_LIBRARIES} pybind11 Threads::Threads)
#target_link_libraries(settestserial ${PYTHON_LIBRARIES} pybind11 Threads::Threads)
#target_link_libraries(settestserialfast ${PYTHON_LIBRARIES} pybind11 Threads::Threads)
#set_target_properties( settestserial PROPERTIES COMPILE_DEFINITIONS "KSET" )
#set_target_properties( settestserialfast PROPERTIES COMPILE_DEFINITIONS "KSET" )
#set_target_properties( settest PROPERTIES COMPILE_DEFINITIONS "KSET" )


# Python library
find_package(PythonLibs ${PYTHON_VERSION} REQUIRED)
find_package(PythonInterp ${PYTHON_VERSION} REQUIRED)
include_directories(${PYTHON_INCLUDE_DIR})


# pybind11 library
add_subdirectory(libs/pybind11-2.2.4)

add_subdirectory(libs/CRoaring)

# set c++14
set (CMAKE_CXX_STANDARD 14)

#add_executable( test_set_driver "${SOURCE_DIR}/TestLosePoint.cc;" ${SOURCES} ${SET_SOURCES} )
#target_link_libraries( test_set_driver ${PYTHON_LIBRARIES} )

#add_executable( test_dict_driver "${SOURCE_DIR}/TestLosePoint.cc;" ${SOURCES} ${DICT_SOURCES} )
#target_link_libraries( test_dict_driver ${PYTHON_LIBRARIES} )

set(PYBIND11_PYTHON_VERSION ${PYTHON_VERSION})
#set(bindings_python_version 2.7)

pybind11_add_module( _Kset ${SOURCES} ${SET_SOURCES} )
set_target_properties( _Kset PROPERTIES COMPILE_DEFINITIONS "KSET" )

pybind11_add_module( _Kdict ${SOURCES} ${DICT_SOURCES} )
set_target_properties( _Kdict PROPERTIES COMPILE_DEFINITIONS "KDICT" )

pybind11_add_module( _Kcounter ${SOURCES} ${COUNTER_SOURCES} )
set_target_properties( _Kcounter PROPERTIES COMPILE_DEFINITIONS "KCOUNTER" )

pybind11_add_module( _Kcolor ${SOURCES} ${COLOR_SOURCES} )
set_target_properties( _Kcolor PROPERTIES COMPILE_DEFINITIONS "KCOLOR" )
target_link_libraries( _Kcolor PRIVATE roaring )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")
        message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")

        if(APPLE)
            #add_custom_command(TARGET test_set_driver POST_BUILD
            #    COMMAND dsymutil test_set_driver
            #    )
            #add_custom_command(TARGET test_dict_driver POST_BUILD
            #    COMMAND dsymutil test_dict_driver
            #    )
            message( "Detected Mac OS with clang, modifying exe with dsymutil for debugging." )
        endif()

        #add_custom_command(TARGET check POST_BUILD
        #    COMMAND dsymutil check
        #    )
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
include_directories(_Kdict "${install_dir}/include")
include_directories(_Kset "${install_dir}/include")
#include_directories(_Kcounter "${install_dir}/include")
#include_directories(_Kcolor "${install_dir}/include")
#include_directories(test_dict_driver "${install_dir}/include")

add_dependencies(_Kdict project_jemalloc)
add_dependencies(_Kset project_jemalloc)
add_dependencies(_Kcounter project_jemalloc)
add_dependencies(_Kcolor project_jemalloc)

add_executable(color_parallel_test "${SOURCE_DIR}/TestColorLoad.cxx;" ${COLOR_SOURCES} ${SOURCES})
set_target_properties(color_parallel_test PROPERTIES COMPILE_DEFINITIONS "KCOLOR")
add_dependencies(color_parallel_test project_jemalloc)
target_link_libraries(color_parallel_test ${PYTHON_LIBRARIES} roaring Threads::Threads pybind11)

#include_directories(settest "${install_dir}/include")
#include_directories(settestserial "${install_dir}/include")
#include_directories(settestserialfast "${install_dir}/include")
#add_dependencies(settest project_jemalloc)
#add_dependencies(settestserial project_jemalloc)
#add_dependencies(settestserialfast project_jemalloc)

