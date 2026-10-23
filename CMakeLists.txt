# The name of our project is "PHOS256".  CMakeLists files in this project can
# refer to the root source directory of the project as ${PHOS256_SOURCE_DIR}
# or as ${CMAKE_SOURCE_DIR} and to the root binary directory of the project as
# ${PHOS256_BINARY_DIR} or ${CMAKE_BINARY_DIR}.

# Check if cmake has the required version
cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

# Set name of our project to "PHOS256". Has to be done
# after check of cmake version
project(PHOS256 CXX)

message(STATUS "CMake platform: ${CMAKE_SYSTEM}")
message(STATUS "Build folder: ${PHOS256_BINARY_DIR}")
message(STATUS "Source folder: ${PHOS256_SOURCE_DIR}")
message(STATUS "Installation folder: ${CMAKE_INSTALL_PREFIX}")

# You can change the build type using
# cmake -DCMAKE_BUILD_TYPE=DEBUG | RELEASE | RELWITHDEBINFO | MINSIZEREL ...
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
endif(NOT CMAKE_BUILD_TYPE)

# Be sure about where libraries and binaries are put
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

# Path to additonal modules
set(CMAKE_MODULE_PATH "${PHOS256_SOURCE_DIR}/cmake")

include(CheckCXXCompilerFlag)

find_package(PYTHIA REQUIRED 8.3.15)
find_package(GEANT4 REQUIRED 11.2.0)
find_package(GEANT4VMC REQUIRED 6.6)
#find_package(CLHEP REQUIRED)


  # ROOT dictionaries and maps
find_package(ROOT REQUIRED)
  include(CMakeROOT)

  # Turn some common warnings into errors
  check_cxx_compiler_flag(-Werror=mismatched-new-delete CXX_COMPILER_HAS_MISMATCHED_NEW_DELETE)
  if(CXX_COMPILER_HAS_MISMATCHED_NEW_DELETE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=mismatched-new-delete")
  endif()
  check_cxx_compiler_flag(-Werror=delete-non-virtual-dtor CXX_COMPILER_DELETE_NON_VIRTUAL_DTOR)
  if(CXX_COMPILER_HAS_MISMATCHED_NEW_DELETE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror=delete-non-virtual-dtor")
  endif()

  # Turn OFF some -Werrors
  check_cxx_compiler_flag(-Wno-error=strict-aliasing CXX_COMPILER_STRICT_ALIASING)
  if(CXX_COMPILER_STRICT_ALIASING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error=strict-aliasing")
  endif()

 # General flags -> Should be moved into a configuration file
  set(CMAKE_POSITION_INDEPENDENT_CODE TRUE)
  # Avoid problems with -fPIE (set automatically by the previous line).
  set(CMAKE_CXX_COMPILE_OPTIONS_PIE "")
  set(CMAKE_C_COMPILE_OPTIONS_PIE "")
  set(CMAKE_Fortran_COMPILE_OPTIONS_PIE "")
  set(CMAKE_NO_SYSTEM_FROM_IMPORTED TRUE)


  add_subdirectory(utils)
  add_subdirectory(dataformats)
  add_subdirectory(simulation)
  add_subdirectory(raw)

