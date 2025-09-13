 ################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             #
 #              GNU Lesser General Public Licence (LGPL) version 3,             #
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
# - Try to find CLHEP
# Once done this will define
#
#  CLHEP_FOUND - system has CLHEP
#  CLHEP_INCLUDE_DIR - the CLHEP include directory
#  CLHEP_LIBRARIES - The libraries needed to use CLHEP
#  CLHEP_DEFINITIONS - Compiler switches required for using CLHEP
#

if (CLHEP_INCLUDE_DIR AND CLHEP_LIBRARY_DIR)
  SET (CLHEP_INCLUDE_DIR CLHEP_INCLUDE_DIR-NOTFOUND)
  SET (CLHEP_LIB_DIR CLHEP_LIB_DIR-NOTFOUND)
  SET (CLHEP_PLISTS_LIB_DIR CLHEP_PLISTS_LIB_DIR-NOTFOUND)
endif (CLHEP_INCLUDE_DIR AND CLHEP_LIBRARY_DIR)

MESSAGE(STATUS "Looking for CLHEP...")


# check if clhep is available in Geant4
If (GEANT4_FOUND)
  SET(GEANT4_LIB ${GEANT4_ROOT}/lib)
  find_path(CLHEP_LIBRARY_DIR NAMES libG4clhep.so GEANT4_LIB
    ${GEANT4_LIBRARY_DIR}
  )
MESSAGE(STATUS "GEANT4_LIBRARY_DIR" ${GEANT4_LIBRARY_DIR} " PATHS" ${${Geant4_DIRS}})
  set(CLHEP_LIBRARIES "-L${CLHEP_LIBRARY_DIR} -lG4clhep")  
MESSAGE(STATUS "CLHEP_LIBRARIES" ${CLHEP_LIBRARIES})

  FIND_PATH(CLHEP_INCLUDE_DIR NAMES CLHEP PATHS
    ${GEANT4_INCLUDE_DIR}
    NO_DEFAULT_PATH
  )

EndIf (GEANT4_FOUND)

#If (CLHEP_INCLUDE_DIR AND CLHEP_LIBRARY_DIR)
If (CLHEP_INCLUDE_DIR OR CLHEP_LIBRARY_DIR)
   set(CLHEP_FOUND TRUE)
EndIf (CLHEP_INCLUDE_DIR OR CLHEP_LIBRARY_DIR)
#EndIf (CLHEP_INCLUDE_DIR AND CLHEP_LIBRARY_DIR)

If (CLHEP_FOUND)
  if (NOT CLHEP_FIND_QUIETLY)
    MESSAGE(STATUS "Looking for CLHEP... - found ${CLHEP_LIBRARY_DIR}")
#    message(STATUS "Found CLHEP: ${CLHEP_LIBRARY_DIR}")
    SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${CLHEP_LIBRARY_DIR})
  endif (NOT CLHEP_FIND_QUIETLY)
Else (CLHEP_FOUND)
  if (CLHEP_FIND_REQUIRED)
    message(FATAL_ERROR "Looking for CLHEP... - Not found")
  endif (CLHEP_FIND_REQUIRED)
EndIf (CLHEP_FOUND)

