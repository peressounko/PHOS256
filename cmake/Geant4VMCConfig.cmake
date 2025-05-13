#------------------------------------------------
# The Geant4 Virtual Monte Carlo package
# Copyright (C) 2014 Ivana Hrivnacova
# All rights reserved.
#
# For the licensing terms see geant4_vmc/LICENSE.
# Contact: root-vmc@cern.ch
#-------------------------------------------------

# Configuration file for CMake build for Geant4 VMC packages.
# It defines the following variables
#  Geant4VMC_PREFIX            - installation prefix
#  Geant4VMC_INCLUDE_DIRS      - include directories for Geant4VMC
#  Geant4VMC_LIBRARIES         - libraries to link against
#  Geant4VMC_USE_G4Root        - build option: with internal G4Root
#  Geant4VMC_USE_EXTERN_G4Root - build option: with external G4Root
#  Geant4VMC_USE_VGM           - build option: with VGM
#  Geant4VMC_USE_GEANT4_UI     - build option: with Geant4 UI drivers
#  Geant4VMC_USE_GEANT4_VIS    - build option: with Geant4 Vis drivers
#  Geant4VMC_USE_GEANT4_G3TOG4 - build option: with Geant4 G3toG4 library
#  G4Root_INCLUDE_DIRS         - include directories for G4Root (if with internal G4Root)
#  G4Root_LIBRARIES            - libraries to link against G4Root (if with internal G4Root)
#  Geant4VMC_CMAKE_INSTALL_LIBDIR - library installation directory
#
# I. Hrivnacova, 13/06/2014

# Compute installation prefix relative to this file
get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_prefix "${_dir}/../.." ABSOLUTE)

# Import options
set(Geant4VMC_PREFIX "${_prefix}")
set(Geant4VMC_USE_G4Root ON)
set(Geant4VMC_USE_EXTERN_G4Root OFF)
set(Geant4VMC_USE_VGM    ON)
set(Geant4VMC_USE_GEANT4_UI  ON)
set(Geant4VMC_USE_GEANT4_VIS ON)
set(Geant4VMC_USE_GEANT4_G3TOG4 OFF)
set(Geant4VMC_CMAKE_INSTALL_LIBDIR lib)

# Import targets
if (Geant4VMC_USE_G4Root AND (NOT Geant4VMC_USE_EXTERN_G4Root))
  # message(STATUS "Importing G4Root targets")
  include("${_prefix}/lib/G4Root-6.6.1/G4RootTargets.cmake")
endif()
include("${_prefix}/lib/Geant4VMC-6.6.1/Geant4VMCTargets.cmake")

# These are IMPORTED targets created by Geant4VMCTargets.cmake
set(Geant4VMC_LIBRARIES geant4vmc)
if(Geant4VMC_USE_G4Root)
  set(Geant4VMC_LIBRARIES "${Geant4VMC_LIBRARIES};g4root")
  set(G4Root_LIBRARIES "${Geant4VMC_LIBRARIES};g4root")
endif()

# Set includes
set(Geant4VMC_INCLUDE_DIRS "${_prefix}/include/geant4vmc")
if(Geant4VMC_USE_G4Root)
  set(Geant4VMC_INCLUDE_DIRS "${Geant4VMC_INCLUDE_DIRS};${_prefix}/include/g4root")
  set(G4Root_INCLUDE_DIRS "${_prefix}/include/g4root")
endif()

# Advanced options
mark_as_advanced(Geant4VMC_USE_EXTERN_G4Root)

