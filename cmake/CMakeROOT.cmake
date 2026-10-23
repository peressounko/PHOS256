# ============================================================================
# ROOT utilities (CMake 3.27+, CMP0144 NEW)
# ============================================================================

include_guard(GLOBAL)

# Проверим, что ROOT найден и переменные определены
if(NOT DEFINED ROOT_CINT OR NOT DEFINED ROOT_LIBDIR)
    message(FATAL_ERROR
        "ROOT utilities: ROOT_CINT / ROOT_LIBDIR не определены. "
        "Убедитесь, что find_package(ROOT) выполнен до включения этого файла.")
endif()

# Нормализуем ROOT_CINT: в новых ROOT это rootcling, в старых — rootcint.
# Убедимся, что это исполняемый файл, а не строка с аргументами.
get_filename_component(_ROOT_CINT_REAL "${ROOT_CINT}" REALPATH)
if(NOT EXISTS "${_ROOT_CINT_REAL}")
    message(FATAL_ERROR "ROOT_CINT='${ROOT_CINT}' не существует.")
endif()

# ============================================================================
# Внутренняя генерация словаря
# ============================================================================
function(_generate_dictionary DNAME LDNAME DHDRS DHDRS_DEPS DINCDIRS)
    # ARGV5 — дополнительные определения (может быть пустым)

    # -I флаги
    set(_include_flags "")
    foreach(_dir IN LISTS DINCDIRS)
        if(_dir)
            list(APPEND _include_flags "-I${_dir}")
        endif()
    endforeach()

    # Определения из директории
    get_directory_property(_dir_defs COMPILE_DEFINITIONS)
    set(_global_defs "")
    foreach(_def IN LISTS _dir_defs)
        string(REPLACE "\"" "\\\"" _def_esc "${_def}")
        list(APPEND _global_defs "-D${_def_esc}")
    endforeach()

    # Дополнительные определения (пятый позиционный аргумент)
    set(_extra_defs "")
    if(ARGC GREATER 5 AND NOT "${ARGV5}" STREQUAL "")
        separate_arguments(_extra_defs UNIX_COMMAND "${ARGV5}")
    endif()

    set(_outputs
        "${CMAKE_CURRENT_BINARY_DIR}/lib${DNAME}.rootmap"
        "${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx"
        "${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}_rdict.pcm"
    )

    # КРИТИЧНО: используем cmake -E env, а не "LD_LIBRARY_PATH=... cmd".
    # Иначе CMake передаёт "LD_LIBRARY_PATH=... rootcling" как один токен
    # и шелл пытается выполнить "-f" как команду.
    add_custom_command(
        OUTPUT  ${_outputs}
        COMMAND ${CMAKE_COMMAND} -E env
                "LD_LIBRARY_PATH=${ROOT_LIBDIR}:$ENV{LD_LIBRARY_PATH}"
                "${_ROOT_CINT_REAL}"
                -f "${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}.cxx"
                -rmf "${CMAKE_CURRENT_BINARY_DIR}/lib${DNAME}.rootmap"
                -rml "lib${DNAME}"
                ${_global_defs}
                ${_extra_defs}
                ${_include_flags}
                ${DHDRS}
                "${LDNAME}"
        DEPENDS ${DHDRS_DEPS} "${LDNAME}" "${_ROOT_CINT_REAL}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        COMMENT "Generating dictionary for ${DNAME}"
        VERBATIM
    )

    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/lib${DNAME}.rootmap"
                  "${CMAKE_CURRENT_BINARY_DIR}/G__${DNAME}_rdict.pcm"
            DESTINATION lib)
endfunction()

# ============================================================================
# Публичная обёртка: разворачивает подкаталоги в списке заголовков
# ============================================================================
function(generate_dictionary DNAME LDNAME DHDRS DINCDIRS)
    set(_dhdrs "")
    set(_daddincdirs "")
    foreach(_itm IN LISTS DHDRS)
        string(FIND "${_itm}" "/" _idx)
        if(_idx GREATER -1)
            get_filename_component(_itmdir  "${_itm}" DIRECTORY)
            get_filename_component(_itmbase "${_itm}" NAME)
            list(APPEND _dhdrs "${_itmbase}")
            list(APPEND _daddincdirs "${CMAKE_CURRENT_SOURCE_DIR}/${_itmdir}")
        else()
            list(APPEND _dhdrs "${_itm}")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES _daddincdirs)

    if(_daddincdirs)
        foreach(_dir IN LISTS _daddincdirs)
            include_directories("${_dir}")
        endforeach()
    endif()

    # Пятый аргумент (EXTRADEFINITIONS) пробрасываем как есть
    _generate_dictionary("${DNAME}" "${LDNAME}"
                         "${_dhdrs}" "${DHDRS}"
                         "${DINCDIRS};${_daddincdirs}"
                         "${ARGV4}")
endfunction()

# ============================================================================
# Генерация .rootmap для ROOT < 6 (rootmap утилита)
# ============================================================================
function(generate_rootmap LIBNAME LIBDEPS LINKDEF)
    if(NOT ROOT_VERSION_MAJOR LESS 6)
        return()
    endif()

    if(NOT DEFINED ROOT_LIBMAP OR NOT EXISTS "${ROOT_LIBMAP}")
        message(WARNING "generate_rootmap: ROOT_LIBMAP не найден, пропуск.")
        return()
    endif()

    set(_local_deps "")
    foreach(_file IN LISTS LIBDEPS)
        get_filename_component(_ext "${_file}" EXT)
        if(_ext)
            list(APPEND _local_deps "${_file}")
        else()
            list(APPEND _local_deps "lib${_file}")
        endif()
    endforeach()

    set(_rootmap "${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap")

    add_custom_command(
        OUTPUT  "${_rootmap}"
        COMMAND ${CMAKE_COMMAND} -E env
                "LD_LIBRARY_PATH=${ROOT_LIBDIR}:$ENV{LD_LIBRARY_PATH}"
                "${ROOT_LIBMAP}"
                -o "${_rootmap}"
                -l "lib${LIBNAME}"
                -d ${_local_deps}
                -c "${LINKDEF}"
        DEPENDS ${LIBNAME} "${LINKDEF}"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        COMMENT "Generating ROOT map for ${LIBNAME}"
        VERBATIM
    )
    add_custom_target(lib${LIBNAME}.rootmap ALL DEPENDS "${_rootmap}")
    install(FILES "${_rootmap}" DESTINATION lib)
endfunction()

# ============================================================================
# Static / helper utilities (без изменений по сути, но аккуратнее)
# ============================================================================
function(generate_static_dependencies shared_list static_list)
    set(_static "")
    foreach(_lib IN LISTS shared_list)
        list(APPEND _static "${_lib}-static")
    endforeach()
    set(${static_list} "${_static}" PARENT_SCOPE)
endfunction()

function(prepend_prefix INLIST PREFIX)
    set(_out "")
    foreach(_item IN LISTS ${INLIST})
        list(APPEND _out "${PREFIX}${_item}")
    endforeach()
    set(${INLIST} "${_out}" PARENT_SCOPE)
endfunction()

function(install_relative)
    set(_mode "FILES")
    set(_files "")
    set(_dest "")
    set(_rest "")
    set(_i 0)
    set(_args ${ARGN})
    list(LENGTH _args _n)

    while(_i LESS _n)
        list(GET _args ${_i} _arg)
        math(EXPR _i "${_i} + 1")

        if(_mode STREQUAL "FILES")
            if(_arg STREQUAL "FILES")
                set(_mode "COLLECT")
            else()
                message(FATAL_ERROR
                    "install_relative() поддерживает только форму install(FILES ...)")
            endif()
        elseif(_mode STREQUAL "COLLECT")
            if(_arg STREQUAL "DESTINATION")
                set(_mode "DEST")
            else()
                list(APPEND _files "${_arg}")
            endif()
        elseif(_mode STREQUAL "DEST")
            set(_dest "${_arg}")
            set(_mode "REST")
        else()
            list(APPEND _rest "${_arg}")
        endif()
    endwhile()

    foreach(_file IN LISTS _files)
        string(SUBSTRING "${_file}" 0 1 _first)
        if(_first STREQUAL "/")
            message(FATAL_ERROR
                "install_relative(): абсолютные пути не поддерживаются (${_file})")
        endif()
        get_filename_component(_prefix "${_file}" DIRECTORY)
        if(_prefix)
            install(FILES "${_file}" DESTINATION "${_dest}/${_prefix}" ${_rest})
        else()
            install(FILES "${_file}" DESTINATION "${_dest}" ${_rest})
        endif()
    endforeach()
endfunction()