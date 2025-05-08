if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.12)
  cmake_policy(SET CMP0074 NEW) # allow to have _ROOT suffix in variables containing program installation roots
endif()

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.27)
  cmake_policy(SET CMP0144 NEW) # search also uppercase variables (terminated by _ROOT) in addition to original case
endif()

include(cmake/CMakeListsDefaults.cmake) # preliminary check of input variables
include(cmake/ROOTMacros.cmake)         # some tools 

enable_language(CXX)

find_package(ROOT REQUIRED)
find_package(Geant4 REQUIRED)
#find_package(GSL REQUIRED)
#find_package(Pythia8 REQUIRED)

#set(BASE_INCLUDE_DIRECTORIES ${BASE_INCLUDE_DIRECTORIES} ${SIMPATH}/include/root ${SIMPATH}/include/vmc)
list(APPEND BASE_INCLUDE_DIR ${ROOT_INCLUDE_DIR} ${Geant4_INCLUDE_DIRS} )

list(APPEND BASE_LIBRARY_DIR ${ROOT_LIBRARY_DIR} )

if(VMC_FOUND)
  message(STATUS "Found VMC: ${VMC_ROOT}")
  list(APPEND BASE_INCLUDE_DIR ${VMC_INCLUDE_DIRS})
  list(APPEND BASE_LIBRARY_DIR ${VMC_LIBDIR})
endif()

#include(cmake/CMakeListsDevelopment.cmake) # developers packages that are not yet compulsory

# START this is for compatibility with legacy, remove later
#set(BASE_INCLUDE_DIRECTORIES ${BASE_INCLUDE_DIR})
#set(BASE_LIBRARY_DIRECTORIES ${BASE_LIBRARY_DIR})
# END this is for compatibility with legacy, remove later

set(CMAKE_INSTALL_LIBDIR "lib")
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# Recurse into the given subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
# LEVEL 1
add_subdirectory (dataformats) # INDEPENDENT
#add_subdirectory (base) # Base
#add_subdirectory (simulation) # INDEPENDENT
#add_subdirectory (reconstruction) # INDEPENDENT


#install(DIRECTORY gconfig/ DESTINATION gconfig PATTERN "gconfig/legacy" EXCLUDE) # copy all except for legacy
#install(DIRECTORY input/ DESTINATION input)
#install(DIRECTORY geometry/ DESTINATION geometry)
#install(DIRECTORY macros/ DESTINATION macros)
#install(FILES scripts/env.sh DESTINATION config)
