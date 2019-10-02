cmake_minimum_required(VERSION 3.1)
include(ExternalProject)

project(run_bft)

add_executable(bft_kmers "run_bft_kmers.cxx")
add_executable(bft_genome "run_bft_genome.cxx")

ExternalProject_Add(project_bft
  GIT_REPOSITORY https://github.com/Colelyman/BloomFilterTrie
  GIT_TAG 0.9.0
  CMAKE_ARGS
  INSTALL_COMMAND "" # doesn't need to be installed system-wide, just locally
)
ExternalProject_Get_Property(project_bft install_dir)
add_library(bft STATIC IMPORTED)
set(BFT_PATH ${install_dir}/src/project_bft)
set_property(TARGET bft PROPERTY IMPORTED_LOCATION ${BFT_PATH}/bin/libbft_static.a)

add_library(judy STATIC IMPORTED)
set(JUDY_PATH ${install_dir}/src/project_bft/bin/judy-1.0.5/lib/libJudy.a)
set_property(TARGET judy PROPERTY IMPORTED_LOCATION ${JUDY_PATH})
#target_link_libraries(bft judy)
include_directories(${BFT_PATH}/include)

add_dependencies(bft_kmers project_bft)
add_dependencies(bft_genome project_bft)
target_link_libraries(bft_kmers bft judy)
target_link_libraries(bft_genome bft judy)
