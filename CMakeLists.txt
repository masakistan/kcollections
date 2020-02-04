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

# BOOST :(
#set(BOOST_REQUESTED_VERSION 1.72.0)
#set(BOOST_ROOT_DIR ${CMAKE_SOURCE_DIR}/libs/boost-${BOOST_REQUESTED_VERSION})
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(Boost ${BOOST_REQUESTED_VERSION} REQUIRED COMPONENTS serialization)
message("boost libs: ${Boost_LIBRARIES} ${Boost_INCLUDE_DIRS}")

# The version number.
set (Kdict_VERSION_MAJOR 1)
set (Kdict_VERSION_MINOR 0)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(SOURCE_DIR "kcollections/src")

message("looking in ${CMAKE_MODULE_PATH}")

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

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
  "${SOURCE_DIR}/globals.cpp"
  "${SOURCE_DIR}/UContainer.cc"
  "${SOURCE_DIR}/Vertex.cc"
  "${SOURCE_DIR}/Kcontainer.cc"
 )

option(PYTHON "compile the python bindings")
if (${PYTHON})
    message("compiling python bindings")
    # pybind11 library
    add_subdirectory(libs/pybind11-2.4.3)
    set(PYBIND11_PYTHON_VERSION ${PYTHON_VERSION})
    include_directories(${PYTHON_INCLUDE_DIR})
    include_directories(libs/pybind11-2.4.3/include)

    pybind11_add_module(_Kset ${kcollections_src} "${SOURCE_DIR}/Kset.cc" "${SOURCE_DIR}/Kcollections.cc" ${uint256_src})
    target_link_libraries(_Kset PRIVATE pybind11::module boost_serialization)
    target_compile_definitions(_Kset PUBLIC "PYTHON" PUBLIC "KSET")

    pybind11_add_module(_Kdict ${kcollections_src} "${SOURCE_DIR}/Kdict.cc" "${SOURCE_DIR}/Kcollections.cc" ${uint256_src})
    target_link_libraries(_Kdict PRIVATE pybind11::module boost_serialization)
    target_compile_definitions(_Kdict PUBLIC "PYTHON" PUBLIC "KDICT")

    pybind11_add_module(_Kcounter ${kcollections_src} "${SOURCE_DIR}/Kcounter.cc" "${SOURCE_DIR}/Kcollections.cc" ${uint256_src})
    target_link_libraries(_Kcounter PRIVATE pybind11::module boost_serialization)
    target_compile_definitions(_Kcounter PUBLIC "PYTHON" PUBLIC "KCOUNTER")
else()
    message("compiling other stuff")
    add_library(kcollections_set_objs OBJECT ${kcollections_src} "${SOURCE_DIR}/Kset.cc")
    target_compile_definitions(kcollections_set_objs PUBLIC "KSET")

    add_library(kcollections_dict_objs OBJECT ${kcollections_src} "${SOURCE_DIR}/Kdict.cc")
    target_compile_definitions(kcollections_dict_objs PUBLIC "KDICT")

    add_library(kcollections_counter_objs OBJECT ${kcollections_src} "${SOURCE_DIR}/Kcounter.cc")
    target_compile_definitions(kcollections_counter_objs PUBLIC "KCOUNTER")

    add_library(Kset STATIC $<TARGET_OBJECTS:kcollections_set_objs> $<TARGET_OBJECTS:uint256_objs>)
	target_link_libraries(Kset boost_serialization)

    add_library(Kdict STATIC $<TARGET_OBJECTS:kcollections_dict_objs> $<TARGET_OBJECTS:uint256_objs>)
	target_link_libraries(Kdict boost_serialization)

    add_library(Kcounter STATIC $<TARGET_OBJECTS:kcollections_counter_objs> $<TARGET_OBJECTS:uint256_objs>)
	target_link_libraries(Kcounter boost_serialization)

    # Testing
    add_executable(counterbulktest "${SOURCE_DIR}/TestCounterBulk.cpp" $<TARGET_OBJECTS:kcollections_counter_objs>)
    target_link_libraries(counterbulktest Threads::Threads Kcounter)
    set_target_properties(counterbulktest PROPERTIES COMPILE_DEFINITIONS "KCOUNTER")

	# Testing kcounter serial
    add_executable(counterserialtest "${SOURCE_DIR}/TestCounterSerial.cpp" $<TARGET_OBJECTS:kcollections_counter_objs>)
    target_link_libraries(counterserialtest Threads::Threads Kcounter)
    set_target_properties(counterserialtest PROPERTIES COMPILE_DEFINITIONS "KCOUNTER")

	# Testing Kdict
    add_executable(dictbulktest "${SOURCE_DIR}/TestDictBulk.cpp" $<TARGET_OBJECTS:kcollections_dict_objs>)
    target_link_libraries(dictbulktest Threads::Threads Kdict)
    set_target_properties(dictbulktest PROPERTIES COMPILE_DEFINITIONS "KDICT")

	# Testing Kset
    add_executable(setbulktest "${SOURCE_DIR}/TestSetBulk.cpp" $<TARGET_OBJECTS:kcollections_set_objs>)
    target_link_libraries(setbulktest Threads::Threads Kset)
    set_target_properties(setbulktest PROPERTIES COMPILE_DEFINITIONS "KSET")
	
	# Testing Kset save
    add_executable(setsavetest "${SOURCE_DIR}/TestSetSave.cpp" $<TARGET_OBJECTS:kcollections_set_objs>)
    target_link_libraries(setsavetest Threads::Threads Kset)
    set_target_properties(setsavetest PROPERTIES COMPILE_DEFINITIONS "KSET")
	
	# Testing Kcounter save
    add_executable(countersavetest "${SOURCE_DIR}/TestCounterSave.cpp" $<TARGET_OBJECTS:kcollections_counter_objs>)
    target_link_libraries(countersavetest Threads::Threads Kset)
    set_target_properties(countersavetest PROPERTIES COMPILE_DEFINITIONS "KCOUNTER")
endif()

# Additional compiler flags
message("compiler id ${CMAKE_CXX_COMPILER_ID}")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    if( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
        string(APPEND CMAKE_CXX_FLAGS_DEBUG " -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")
        message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
        if(APPLE)
            add_custom_command(TARGET counterbulktest POST_BUILD
                COMMAND dsymutil dictbulktest
                )
             add_custom_command(TARGET counterbulktest POST_BUILD
                COMMAND dsymutil counterbulktest
                )
            message( "Detected Mac OS with clang, modifying exe with dsymutil for debugging." )
        endif()
    endif()
  # using Clang
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  #string( APPEND CMAKE_CXX_FLAGS " -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free" )
  message( "Using CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}" )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
endif()
