# The name of our project is "PHOS256".  CMakeLists files in this project can
# refer to the root source directory of the project as ${PHOS256_SOURCE_DIR}
# or as ${CMAKE_SOURCE_DIR} and to the root binary directory of the project as
# ${PHOS256_BINARY_DIR} or ${CMAKE_BINARY_DIR}.

# Check if cmake has the required version
cmake_minimum_required(VERSION 3.20.0 FATAL_ERROR)

# Set name of our project to "PHOS256". Has to be done
# after check of cmake version
project(PHOS256 NONE)

include(cmake/CMakeLists.cmake)
