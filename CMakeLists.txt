cmake_minimum_required (VERSION 3.1)
include(ExternalProject)

# Fix behavior of CMAKE_CXX_STANDARD when targeting macOS.
if (POLICY CMP0025)
    cmake_policy(SET CMP0025 NEW)
endif ()

set(CMAKE_VERBOSE_MAKEFILE on)
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
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(SOURCE_DIR "kcollections/src")

message("looking in ${CMAKE_MODULE_PATH}")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(CMAKE_CXX_STANDARD 14)

# Add inc directories
include_directories(inc)

include_directories(libs/uint256_t)

# Add src files
set(uint256_src
  "libs/uint256_t/uint128_t.cpp"
  "libs/uint256_t/uint256_t.cpp"
  )
add_library(uint256_objs OBJECT ${uint256_src})

set(kcollections_src
  "${SOURCE_DIR}/MurmurHash3.cc"
  "${SOURCE_DIR}/UContainer.cc"
  "${SOURCE_DIR}/Vertex.cc"
  "${SOURCE_DIR}/Kcontainer.cc"
 )

option(PYTHON "compile the python bindings")
if (${PYTHON})
    # pybind11 library
    add_subdirectory(libs/pybind11-2.3.0)
    set(PYBIND11_PYTHON_VERSION ${PYTHON_VERSION})
    include_directories(${PYTHON_INCLUDE_DIR})
    include_directories(libs/pybind11-2.3.0/include)

    pybind11_add_module(_Kset ${kcollections_src} "${SOURCE_DIR}/Kset.cc" "${SOURCE_DIR}/Kcollections.cc" ${uint256_src})
    target_link_libraries(_Kset PRIVATE pybind11::module)
    target_compile_definitions(_Kset PUBLIC "PYTHON" PUBLIC "KSET")
    # add_dependencies(_Kset project_jemalloc)

    pybind11_add_module(_Kdict ${kcollections_src} "${SOURCE_DIR}/Kdict.cc" "${SOURCE_DIR}/Kcollections.cc" ${uint256_src})#$<TARGET_OBJECTS:python_cxx_objs>)
    target_link_libraries(_Kdict PRIVATE pybind11::module)
    target_compile_definitions(_Kdict PUBLIC "PYTHON" PUBLIC "KDICT")
    # add_dependencies(_Kdict project_jemalloc)

    pybind11_add_module(_Kcounter ${kcollections_src} "${SOURCE_DIR}/Kcounter.cc" "${SOURCE_DIR}/Kcollections.cc" ${uint256_src})#$<TARGET_OBJECTS:kcollections_counter_objs> $<TARGET_OBJECTS:python_cxx_objs>)
    target_link_libraries(_Kcounter PRIVATE pybind11::module)
    target_compile_definitions(_Kcounter PUBLIC "PYTHON" PUBLIC "KCOUNTER")
    # add_dependencies(_Kcounter project_jemalloc)
else()
    add_library(kcollections_set_objs OBJECT ${kcollections_src} "${SOURCE_DIR}/Kset.cc")
    target_compile_definitions(kcollections_set_objs PUBLIC "KSET")

    add_library(kcollections_dict_objs OBJECT ${kcollections_src} "${SOURCE_DIR}/Kdict.cc")
    target_compile_definitions(kcollections_dict_objs PUBLIC "KDICT")

    add_library(kcollections_counter_objs OBJECT ${kcollections_src} "${SOURCE_DIR}/Kcounter.cc")
    target_compile_definitions(kcollections_counter_objs PUBLIC "KCOUNTER")

    add_library(Kset STATIC $<TARGET_OBJECTS:kcollections_set_objs> $<TARGET_OBJECTS:uint256_objs>)
    # add_dependencies(Kset project_jemalloc)
    # target_link_libraries(Kset jemalloc)

    add_library(Kdict STATIC $<TARGET_OBJECTS:kcollections_dict_objs> $<TARGET_OBJECTS:uint256_objs>)
    # add_dependencies(Kdict project_jemalloc)
    # target_link_libraries(Kdict jemalloc)

    add_library(Kcounter STATIC $<TARGET_OBJECTS:kcollections_counter_objs> $<TARGET_OBJECTS:uint256_objs>)
    # add_dependencies(Kcounter project_jemalloc)
    # target_link_libraries(Kcounter jemalloc)

    # Testing
    add_executable(counterbulktest "${SOURCE_DIR}/TestCounterBulk.cpp" $<TARGET_OBJECTS:kcollections_counter_objs>)
    target_link_libraries(counterbulktest Threads::Threads Kcounter)
    set_target_properties(counterbulktest PROPERTIES COMPILE_DEFINITIONS "KCOUNTER")
    # include_directories(counterbulktest "${install_dir}/include")
    # add_dependencies(counterbulktest project_jemalloc)
endif()

# Additinoal compiler flags
message("compiler id ${CMAKE_CXX_COMPILER_ID}")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")
        message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
        if(APPLE)
            add_custom_command(TARGET counterbulktest POST_BUILD
                COMMAND dsymutil counterbulktest
                )
            message( "Detected Mac OS with clang, modifying exe with dsymutil for debugging." )
        endif()
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
# set(PROJECT_JEMALLOC_PATH ${CMAKE_BINARY_DIR}/jemalloc-5.2.2)
# ExternalProject_Add(project_jemalloc
#   URL https://github.com/jemalloc/jemalloc/releases/download/5.2.1/jemalloc-5.2.1.tar.bz2
#   URL_HASH MD5=3d41fbf006e6ebffd489bdb304d009ae
#   PREFIX ${PROJECT_JEMALLOC_PATH}
#   CONFIGURE_COMMAND ${PROJECT_JEMALLOC_PATH}/src/project_jemalloc/configure --prefix=${PROJECT_JEMALLOC_PATH}
#   STEP_TARGETS build
#   EXCLUDE_FROM_ALL TRUE
#   BUILD_IN_SOURCE TRUE
#   BUILD_COMMAND ${MAKE})
# ExternalProject_Get_Property(project_jemalloc install_dir)
# add_library(jemalloc STATIC IMPORTED)
# add_library(jemalloc_shared SHARED IMPORTED)
# set_property(TARGET jemalloc PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libjemalloc${CMAKE_STATIC_LIBRARY_SUFFIX})
# set_property(TARGET jemalloc PROPERTY INTERFACE_LINK_LIBRARIES pthread dl)
# set_property(TARGET jemalloc_shared PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libjemalloc${CMAKE_SHARED_LIBRARY_SUFFIX})

# include_directories(Kset "${install_dir}/include")
# include_directories(Kdict "${install_dir}/include")
# include_directories(Kcounter "${install_dir}/include")

# if (${PYTHON})
#   include_directories(_Kset "${install_dir}/include")
#   include_directories(_Kdict "${install_dir}/include")
#   include_directories(_Kcounter "${install_dir}/include")
# endif()
