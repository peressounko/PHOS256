set(PTH_FOUND "")
set(PTH_MISSING "")
set(PTH_ADDEDNDUM "")

macro(checkPath defaultName envName description)
  if (NOT "${${defaultName}}" STREQUAL "") # I have already set variable with path
    list(APPEND PTH_FOUND "${description}: ${${defaultName}}")
  elseif (NOT "$ENV{${envName}}" STREQUAL "") # I have global variable of expected name
    set(${defaultName} $ENV{${envName}})
    list(APPEND PTH_FOUND "${description}: ${${defaultName}}")
  else () # I don't have global variable of given name
    list(APPEND PTH_MISSING "${description}: -D${defaultName} (variable ${envName})")
    string(APPEND PTH_ADDENDUM "\nexport ${defaultName}=<path to ${envName}>")
  endif()
endmacro()

#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()


#if ("$ENV{PHOS256_ROOT}" STREQUAL "") # user did not set PHOS256_ROOT variable
#  message(FATAL_ERROR "$PHOS256_ROOT not set.\n Build without setting environment variable holding installation (and execution) directory is not allowed.")
#endif()

# force CMAKE_INSTALL_PREFIX to PHOS256_ROOT
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX $ENV{PHOS256_ROOT} CACHE PATH "" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

checkPath(CMAKE_INSTALL_PREFIX PHOS256_ROOT "PHOS256_ROOT")
#checkPath(BOOST_ROOT BOOST_ROOT "Boost root")
checkPath(GEANT4_ROOT GEANT4_ROOT "Geant4 root")

#if ((NOT ("${PC_LIBXML_INCLUDEDIR}" STREQUAL "")) AND (NOT ("${PC_LIBXML_LIBDIR}" STREQUAL ""))) # if user has provided LIBXML path, don't look for LIBXML2_ROOT
#  list(APPEND PTH_FOUND "libXML2 include dir: ${PC_LIBXML_INCLUDEDIR}")
#  list(APPEND PTH_FOUND "libXML2 library dir: ${PC_LIBXML_LIBDIR}")
#else()
#  checkPath(LIBXML2_ROOT LIBXML2_ROOT "libXML2 root")
#endif()
#checkPath(FMT_ROOT FMT_ROOT "FMT root")

if (PTH_MISSING) # Some dependencies are missing. I will write to the screen, what I have found and what I am missing
  if (PTH_FOUND)
    message("${Green}\nPaths (packages) found:${ColourReset}")
    foreach (_PTH_FOUND ${PTH_FOUND})
      message("${BoldGreen}${_PTH_FOUND}${ColourReset}")
    endforeach()
  endif()

  if (PTH_MISSING)
    message("${Red}\nPaths (packages) missing:${ColourReset}")
    foreach (_PTH_MISSING ${PTH_MISSING})
      message("${BoldRed}${_PTH_MISSING}${ColourReset}")
    endforeach()
    message("${Yellow}\nCreate the missing variables:${ColourReset}")
    message("${BoldYellow}${PTH_ADDENDUM}\n${ColourReset}")
    message("${Yellow}and re-run your cmake command${ColourReset}\n")

    # go through input commandline parameters and build proper commandline string. Not necessary, just polite.
    get_cmake_property(CACHE_VARS CACHE_VARIABLES)
    foreach(CACHE_VAR ${CACHE_VARS})
      get_property(CACHE_VAR_HELPSTRING CACHE ${CACHE_VAR} PROPERTY HELPSTRING)
      if(CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line.")
        get_property(CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE)
        if(CACHE_VAR_TYPE STREQUAL "UNINITIALIZED")
          set(CACHE_VAR_TYPE)
        else()
          set(CACHE_VAR_TYPE :${CACHE_VAR_TYPE})
        endif()
        set(CMAKE_ARGS "${CMAKE_ARGS} -D${CACHE_VAR}${CACHE_VAR_TYPE}=\"${${CACHE_VAR}}\"")
      endif()
    endforeach()
    message("${BoldYellow}cmake${CMAKE_ARGS} ${CMAKE_SOURCE_DIR}${ColourReset}\n")
    message(FATAL_ERROR "") # terminate function
  endif()
else() # I have found all dependencies, no text to the screen
  if ("${ENV${Geant3_DIR}}" STREQUAL "") # I have not set Geant3_DIR
    set(ENV{Geant3_DIR} ${GEANT3_ROOT})
  endif()
  set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${LIBXML2_ROOT}/lib/pkgconfig")
endif()

if ("${VMC_ROOT}" STREQUAL "") # if VMC_ROOT was not set, copy it from ENV if possible
  if (NOT "$ENV{VMC_ROOT}" STREQUAL "")
    set(VMC_ROOT $ENV{VMC_ROOT})
  endif()
endif()

    # Check for rootcint
    find_program(ROOT_CINT NAMES rootcint PATHS ${ROOTSYS}/bin NO_DEFAULT_PATH)
    if(ROOT_CINT)
        message(STATUS "Found ${ROOT_CINT}")
    else()
        message(FATAL_ERROR "Could not find rootcint executable.")
    endif(ROOT_CINT)


if (CMAKE_INSTALL_PREFIX STREQUAL CMAKE_SOURCE_DIR) # user wants to install PHOS256_ROOT into sources directory
  message(FATAL_ERROR "${BoldRed}\n$PHOS256_ROOT may not point to the source directory.${ColourReset}\n${BoldYellow}In-source install is not allowed.${ColourReset}\n")
endif()
