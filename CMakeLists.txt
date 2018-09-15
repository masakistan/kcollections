cmake_minimum_required (VERSION 2.6)
project (Kdict)
# The version number.
set (Kdict_VERSION_MAJOR 0)
set (Kdict_VERSION_MINOR 1)

# Add inc directory
include_directories(inc)

# Add src files
file(GLOB SOURCES "src/*.cc")

# Python library
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# set c++11
set (CMAKE_CXX_STANDARD 11)

# configure a header file to pass some of the CMake settings
# to the source code
#configure_file (
#  "${PROJECT_SOURCE_DIR}/KdictConfig.h.in"
#  "${PROJECT_BINARY_DIR}/KdictConfig.h"
#  )

add_executable(Kdict ${SOURCES})
target_link_libraries(Kdict ${PYTHON_LIBRARIES})
