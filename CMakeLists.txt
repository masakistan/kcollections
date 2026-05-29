cmake_minimum_required(VERSION 3.18)
project(kcollections LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# pybind11: prefer the copy installed by the build backend
find_package(pybind11 CONFIG REQUIRED)

find_package(Threads REQUIRED)

include_directories(${CMAKE_SOURCE_DIR}/inc)
include_directories(${CMAKE_SOURCE_DIR}/libs/uint256_t)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SOURCE_DIR "${CMAKE_SOURCE_DIR}/kcollections/src")

set(uint256_src
  "${CMAKE_SOURCE_DIR}/libs/uint256_t/uint128_t.cpp"
  "${CMAKE_SOURCE_DIR}/libs/uint256_t/uint256_t.cpp"
)

set(kcollections_src
  "${SOURCE_DIR}/globals.cpp"
  "${SOURCE_DIR}/UContainer.cc"
  "${SOURCE_DIR}/Vertex.cc"
  "${SOURCE_DIR}/Kcontainer.cc"
)

function(add_kcollections_module target_name main_src compile_def)
  pybind11_add_module(${target_name} MODULE
    ${kcollections_src}
    ${main_src}
    "${SOURCE_DIR}/Kcollections.cc"
    ${uint256_src}
  )
  target_compile_definitions(${target_name} PRIVATE PYTHON ${compile_def})
  target_link_libraries(${target_name} PRIVATE
    pybind11::module
    Threads::Threads
  )
  target_include_directories(${target_name} PRIVATE
    ${CMAKE_SOURCE_DIR}/inc
    ${CMAKE_SOURCE_DIR}/libs/uint256_t
  )
  install(TARGETS ${target_name} LIBRARY DESTINATION kcollections)
endfunction()

add_kcollections_module(_Kset "${SOURCE_DIR}/Kset.cc" KSET)
add_kcollections_module(_Kdict "${SOURCE_DIR}/Kdict.cc" KDICT)
add_kcollections_module(_Kcounter "${SOURCE_DIR}/Kcounter.cc" KCOUNTER)

option(KCOLLECTIONS_BUILD_CPP_TESTS "Build C++ unit test binaries" OFF)
if(KCOLLECTIONS_BUILD_CPP_TESTS)
  add_library(uint256_objs OBJECT ${uint256_src})

  macro(add_cpp_test exe_name extra_src compile_def link_lib)
    add_executable(${exe_name}
      ${kcollections_src}
      ${extra_src}
      $<TARGET_OBJECTS:uint256_objs>
    )
    target_compile_definitions(${exe_name} PRIVATE ${compile_def})
    target_include_directories(${exe_name} PRIVATE
      ${CMAKE_SOURCE_DIR}/inc
      ${CMAKE_SOURCE_DIR}/libs/uint256_t
    )
    target_link_libraries(${exe_name} PRIVATE Threads::Threads ${link_lib})
  endmacro()

  add_library(kcollections_set_objs OBJECT ${kcollections_src} "${SOURCE_DIR}/Kset.cc")
  target_compile_definitions(kcollections_set_objs PRIVATE KSET)
  add_library(kcollections_dict_objs OBJECT ${kcollections_src} "${SOURCE_DIR}/Kdict.cc")
  target_compile_definitions(kcollections_dict_objs PRIVATE KDICT)
  add_library(kcollections_counter_objs OBJECT ${kcollections_src} "${SOURCE_DIR}/Kcounter.cc")
  target_compile_definitions(kcollections_counter_objs PRIVATE KCOUNTER)

  add_library(Kset STATIC $<TARGET_OBJECTS:kcollections_set_objs> $<TARGET_OBJECTS:uint256_objs>)
  add_library(Kdict STATIC $<TARGET_OBJECTS:kcollections_dict_objs> $<TARGET_OBJECTS:uint256_objs>)
  add_library(Kcounter STATIC $<TARGET_OBJECTS:kcollections_counter_objs> $<TARGET_OBJECTS:uint256_objs>)

  add_executable(setbulktest "${SOURCE_DIR}/TestSetBulk.cpp")
  target_link_libraries(setbulktest Threads::Threads Kset)
  target_compile_definitions(setbulktest PRIVATE KSET)

  add_executable(dictbulktest "${SOURCE_DIR}/TestDictBulk.cpp")
  target_link_libraries(dictbulktest Threads::Threads Kdict)
  target_compile_definitions(dictbulktest PRIVATE KDICT)

  add_executable(counterbulktest "${SOURCE_DIR}/TestCounterBulk.cpp")
  target_link_libraries(counterbulktest Threads::Threads Kcounter)
  target_compile_definitions(counterbulktest PRIVATE KCOUNTER)

  add_executable(setsavetest "${SOURCE_DIR}/TestSetSave.cpp")
  target_link_libraries(setsavetest Threads::Threads Kset)
  target_compile_definitions(setsavetest PRIVATE KSET)

  add_executable(dictsavetest "${SOURCE_DIR}/TestDictSave.cpp")
  target_link_libraries(dictsavetest Threads::Threads Kdict)
  target_compile_definitions(dictsavetest PRIVATE KDICT)

  add_executable(countersavetest "${SOURCE_DIR}/TestCounterSave.cpp")
  target_link_libraries(countersavetest Threads::Threads Kcounter)
  target_compile_definitions(countersavetest PRIVATE KCOUNTER)
endif()
