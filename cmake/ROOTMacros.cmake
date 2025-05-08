 ################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             #
 #              GNU Lesser General Public Licence (LGPL) version 3,             #
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
Function(Format _output input prefix suffix)

# DevNotes - input should be put in quotes or the complete list does not get passed to the function
  set(format)
  foreach(arg ${input})
    set(item ${arg})
    if(prefix)
      string(REGEX MATCH "^${prefix}" pre ${arg})
    endif(prefix)
    if(suffix)
      string(REGEX MATCH "${suffix}$" suf ${arg})
    endif(suffix)
    if(NOT pre)
      set(item "${prefix}${item}")
    endif(NOT pre)
    if(NOT suf)
      set(item "${item}${suffix}")
    endif(NOT suf)
    list(APPEND format ${item})
  endforeach(arg)
  set(${_output} ${format} PARENT_SCOPE)

endfunction(Format)


# Generation of the dictionaries
# @DNAME  Dictionary name
# @LDNAME LinkDef file name, ex: LinkDef.h
# @DHDRS  Dictionary headers
# @DHDRS_DEPS  Dictionary header files used as dependencies to the rootmap target
# @DINCDIR Include folders that need to be passed to cint/cling
# @EXTRADEFINITIONS - optional, extra compile flags specific to library
#       - used as ${ARGV5}
macro(_generate_dictionary DNAME LDNAME DHDRS DHDRS_DEPS DINCDIRS)

    # Creating the INCLUDE path for cint/cling
    foreach(_dir ${DINCDIRS})
        set(INCLUDE_PATH -I${_dir} ${INCLUDE_PATH})
    endforeach()

set(DESTINATION "/home/prsnko/tmp/test")
    # Get the list of definitions from the directory to be sent to CINT
    get_directory_property(tmpdirdefs COMPILE_DEFINITIONS)
    foreach(dirdef ${tmpdirdefs})
        string(REPLACE "\"" "\\\"" dirdef_esc ${dirdef})
        set(GLOBALDEFINITIONS -D${dirdef_esc} ${GLOBALDEFINITIONS})
    endforeach()

    # Custom definitions specific to library
    # Received as the forth optional argument
    separate_arguments(EXTRADEFINITIONS UNIX_COMMAND "${ARGV5}")

    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib${DNAME}.rootmap ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}_rdict.pcm
                     COMMAND
                       LD_LIBRARY_PATH=${ROOT_LIBDIR}:$ENV{LD_LIBRARY_PATH} ${ROOT_CINT}
                     ARGS
                       -f ${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx
                       -rmf ${CMAKE_CURRENT_BINARY_DIR}/lib${DNAME}.rootmap -rml lib${DNAME}
                       ${GLOBALDEFINITIONS} ${EXTRADEFINITIONS} ${INCLUDE_PATH} ${DHDRS} ${LDNAME}
                     DEPENDS
                       ${DHDRS_DEPS} ${LDNAME} ${ROOT_CINT}
                     WORKING_DIRECTORY
                       ${CMAKE_CURRENT_BINARY_DIR}
                    )

  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/lib${DNAME}.rootmap" DESTINATION lib)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}_rdict.pcm" DESTINATION lib)

endmacro(_generate_dictionary)


# Same as generate_dictionary, but flattens the list of headers and sets additional include paths
# with include_directories
macro(generate_dictionary DNAME LDNAME DHDRS DINCDIRS)

    set(_dhdrs "")
    set(_daddincdirs "")
    foreach(_itm ${DHDRS})
      string(FIND "${_itm}" "/" _idx)
      if(_idx GREATER -1)
        # Has a subdirectory specified
        get_filename_component(_itmdir "${_itm}" DIRECTORY)
        get_filename_component(_itmbase "${_itm}" NAME)
        list(APPEND _dhdrs "${_itmbase}")
        list(APPEND _daddincdirs "${CMAKE_CURRENT_SOURCE_DIR}/${_itmdir}")
      else()
        # No subdirectory specified
        list(APPEND _dhdrs "${_itm}")
      endif()
    endforeach()
    list(REMOVE_DUPLICATES _daddincdirs)
    if(NOT "${_daddincdirs}" STREQUAL "")
      foreach(_dir "${_daddincdirs}")
        include_directories("${_dir}")
      endforeach()
    endif()

    _generate_dictionary("${DNAME}" "${LDNAME}" "${_dhdrs}" "${DHDRS}" "${DINCDIRS};${_daddincdirs}" "${ARGV4}")

endmacro(generate_dictionary)


# Generate the ROOTmap files
# @LIBNAME - library name: libAnalysis.so -> Analysis.rootmap
# @LIBDEPS - library dependencies
# @LINKDEF - LinkDef header
macro(generate_rootmap LIBNAME LIBDEPS LINKDEF)
#    message(STATUS "LIBNAME = ${LIBNAME}")
#    message(STATUS "LIBDEPS = ${LIBDEPS}")
#    message(STATUS "LINKDEF = ${LINKDEF}")
#    message(STATUS "ROOT_LIBMAP=${ROOT_LIBMAP}")

if (ROOT_VERSION_MAJOR LESS 6)

    set(LOCAL_DEPS)
    foreach(file ${LIBDEPS})
        get_filename_component(ext ${file} EXT)
        if(ext)
            set(LOCAL_DEPS ${LOCAL_DEPS} ${file})
        else()
            set(LOCAL_DEPS ${LOCAL_DEPS} lib${file})
        endif()
    endforeach()

#    message(STATUS "Generating ROOT map for ${LIBNAME}")
    add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap
                       COMMAND LD_LIBRARY_PATH=${ROOT_LIBDIR}:$ENV{LD_LIBRARY_PATH} ${ROOT_LIBMAP}
                       ARGS -o ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap -l lib${LIBNAME} -d ${LOCAL_DEPS} -c ${LINKDEF}
                       DEPENDS ${LIBNAME}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} VERBATIM
                      )
    add_custom_target(lib${LIBNAME}.rootmap ALL DEPENDS  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap DESTINATION lib)

endif (ROOT_VERSION_MAJOR LESS 6)

endmacro(generate_rootmap)

MACRO (GENERATE_ROOT_TEST_SCRIPT SCRIPT_FULL_NAME)

  get_filename_component(path_name ${SCRIPT_FULL_NAME} PATH)
  get_filename_component(file_extension ${SCRIPT_FULL_NAME} EXT)
  get_filename_component(file_name ${SCRIPT_FULL_NAME} NAME_WE)
  set(shell_script_name "${file_name}.sh")

  #MESSAGE("PATH: ${path_name}")
  #MESSAGE("Ext: ${file_extension}")
  #MESSAGE("Name: ${file_name}")
  #MESSAGE("Shell Name: ${shell_script_name}")

  string(REPLACE ${PROJECT_SOURCE_DIR}
         ${PROJECT_BINARY_DIR} new_path ${path_name}
        )

  #MESSAGE("New PATH: ${new_path}")

  file(MAKE_DIRECTORY ${new_path}/data)

  CONVERT_LIST_TO_STRING(${LD_LIBRARY_PATH})
  set(MY_LD_LIBRARY_PATH ${output})

  CONVERT_LIST_TO_STRING(${ROOT_INCLUDE_PATH})
  set(MY_ROOT_INCLUDE_PATH ${output})

  set(my_script_name ${SCRIPT_FULL_NAME})

  Write_Geant4Data_Variables_sh()
  IF(FAIRROOTPATH)
    configure_file(${FAIRROOTPATH}/share/fairbase/cmake/scripts/root_macro.sh.in
                   ${new_path}/${shell_script_name}
                  )
  ELSE(FAIRROOTPATH)
    configure_file(${PROJECT_SOURCE_DIR}/cmake/scripts/root_macro.sh.in
                   ${new_path}/${shell_script_name}
                  )
  ENDIF(FAIRROOTPATH)
  execute_process(COMMAND /bin/chmod u+x ${new_path}/${shell_script_name} OUTPUT_QUIET)

ENDMACRO (GENERATE_ROOT_TEST_SCRIPT)


Macro(ROOT_GENERATE_ROOTMAP)

  # All Arguments needed for this new version of the macro are defined
  # in the parent scope, namely in the CMakeLists.txt of the submodule
  if (DEFINED LINKDEF)
    foreach(l ${LINKDEF})
      If( IS_ABSOLUTE ${l})
        Set(Int_LINKDEF ${Int_LINKDEF} ${l})
      Else( IS_ABSOLUTE ${l})
        Set(Int_LINKDEF ${Int_LINKDEF} ${CMAKE_CURRENT_SOURCE_DIR}/${l})
      EndIf( IS_ABSOLUTE ${l})
    endforeach()

    foreach(d ${DEPENDENCIES})
      get_filename_component(_ext ${d} EXT)
      If(NOT _ext MATCHES a$)
        if(_ext)
          set(Int_DEPENDENCIES ${Int_DEPENDENCIES} ${d})
        else()
          set(Int_DEPENDENCIES ${Int_DEPENDENCIES} lib${d}.so)
        endif()
      Else()
        Message("Found Static library with extension ${_ext}")
      EndIf()
    endforeach()

    set(Int_LIB ${LIBRARY_NAME})
    set(Int_OUTFILE ${LIBRARY_OUTPUT_PATH}/lib${Int_LIB}.rootmap)

    add_custom_command(OUTPUT ${Int_OUTFILE}
                       COMMAND ${RLIBMAP_EXECUTABLE} -o ${Int_OUTFILE} -l ${Int_LIB}
                               -d ${Int_DEPENDENCIES} -c ${Int_LINKDEF}
                       DEPENDS ${Int_LINKDEF} ${RLIBMAP_EXECUTABLE} )
    add_custom_target( lib${Int_LIB}.rootmap ALL DEPENDS  ${Int_OUTFILE})
    set_target_properties(lib${Int_LIB}.rootmap PROPERTIES FOLDER RootMaps )
    #---Install the rootmap file------------------------------------
    #install(FILES ${Int_OUTFILE} DESTINATION lib COMPONENT libraries)
    install(FILES ${Int_OUTFILE} DESTINATION lib)
  endif(DEFINED LINKDEF)
EndMacro(ROOT_GENERATE_ROOTMAP)

Macro(GENERATE_LIBRARY)

  set(Int_LIB ${LIBRARY_NAME})
  Set(HeaderRuleName "${Int_LIB}_HEADER_RULES")
  Set(DictName "G__${Int_LIB}Dict.cxx")

  If(NOT DICTIONARY)
    Set(DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/${DictName})
  EndIf(NOT DICTIONARY)

  If( IS_ABSOLUTE ${DICTIONARY})
    Set(DICTIONARY ${DICTIONARY})
  Else( IS_ABSOLUTE ${DICTIONARY})
    Set(Int_DICTIONARY ${CMAKE_CURRENT_SOURCE_DIR}/${DICTIONARY})
  EndIf( IS_ABSOLUTE ${DICTIONARY})

  Set(Int_SRCS ${SRCS})

  If(HEADERS)
    Set(HDRS ${HEADERS})
  Else(HEADERS)
    CHANGE_FILE_EXTENSION(*.cxx *.h HDRS "${SRCS}")
  EndIf(HEADERS)

  If(IWYU_FOUND)
    Set(_INCLUDE_DIRS ${INCLUDE_DIRECTORIES} ${SYSTEM_INCLUDE_DIRECTORIES})
    CHECK_HEADERS("${Int_SRCS}" "${_INCLUDE_DIRS}" ${HeaderRuleName})
  EndIf(IWYU_FOUND)

  #  install(FILES ${HDRS} DESTINATION include)
  # fix to keep structure of header files
  foreach ( file ${HDRS} )
    cmake_path(GET file PARENT_PATH file_dir)
    cmake_path(GET file FILENAME file_name)
    if("${file_dir}" STREQUAL "") # I have only header file - no (sub)directory
      install(FILES ${file} DESTINATION include)
    else()
      install(FILES ${file} DESTINATION include) # I copy header twice
      install(FILES ${file} DESTINATION include/${file_dir})
    endif()
  endforeach()

 If(LINKDEF)
    If( IS_ABSOLUTE ${LINKDEF})
      Set(Int_LINKDEF ${LINKDEF})
    Else( IS_ABSOLUTE ${LINKDEF})
      Set(Int_LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/${LINKDEF})
    EndIf( IS_ABSOLUTE ${LINKDEF})
    generate_dictionary("${Int_LIB}" "${DictName}" "${HEADERS}" "${_INCLUDE_DIRS}")
    SET(Int_SRCS ${Int_SRCS} ${DICTIONARY})
    SET_SOURCE_FILES_PROPERTIES(${DICTIONARY}
       PROPERTIES COMPILE_FLAGS "-Wno-old-style-cast"
    )
  EndIf(LINKDEF)

  set(Int_DEPENDENCIES)
  foreach(d ${DEPENDENCIES})
    get_filename_component(_ext ${d} EXT)
    If(NOT _ext MATCHES a$)
      set(Int_DEPENDENCIES ${Int_DEPENDENCIES} ${d})
    Else()
      Message("Found Static library with extension ${_ext}")
      get_filename_component(_lib ${d} NAME_WE)
      set(Int_DEPENDENCIES ${Int_DEPENDENCIES} ${_lib})
    EndIf()
  endforeach()

  ############### build the library #####################
  If(${CMAKE_GENERATOR} MATCHES Xcode)
    Add_Library(${Int_LIB} SHARED ${Int_SRCS} ${NO_DICT_SRCS} ${HDRS} ${LINKDEF})
  Else()
    Add_Library(${Int_LIB} SHARED ${Int_SRCS} ${NO_DICT_SRCS} ${LINKDEF})
  EndIf()
  target_link_libraries(${Int_LIB} ${Int_DEPENDENCIES})
#  set_target_properties(${Int_LIB} PROPERTIES ${PROJECT_LIBRARY_PROPERTIES})

  ############### install the library ###################
  install(TARGETS ${Int_LIB} DESTINATION lib)

  Set(LIBRARY_NAME)
  Set(DICTIONARY)
  Set(LINKDEF)
  Set(SRCS)
  Set(HEADERS)
  Set(NO_DICT_SRCS)
  Set(DEPENDENCIES)
EndMacro(GENERATE_LIBRARY)


Macro(GENERATE_EXECUTABLE)

#  If(IWYU_FOUND)
#    Set(HeaderRuleName "${EXE_NAME}_HEADER_RULES")
#    CHECK_HEADERS("${SRCS}" "${INCLUDE_DIRECTORIES}" ${HeaderRuleName})
#  EndIf(IWYU_FOUND)

  ############### build the library #####################
  Add_Executable(${EXE_NAME} ${SRCS})
  target_link_libraries(${EXE_NAME} ${DEPENDENCIES})

  ############### install the library ###################
  if(DEFINED BIN_DESTINATION)
    install(TARGETS ${EXE_NAME} DESTINATION ${BIN_DESTINATION})
  else(DEFINED BIN_DESTINATION)
    install(TARGETS ${EXE_NAME} DESTINATION bin)
  endif(DEFINED BIN_DESTINATION)

  Set(EXE_NAME)
  Set(SRCS)
  Set(DEPENDENCIES)

EndMacro(GENERATE_EXECUTABLE)
