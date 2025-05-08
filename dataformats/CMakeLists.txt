# Create a library called "libDataFormats" which includes the source files given in the array .
Set(LinName "DataFormats")

set(INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_LIST_DIR}
)

Set(SYSTEM_INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIR}
)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
  ${ROOT_LIBRARY_DIR}
)

link_directories(${LINK_DIRECTORIES})

# List of source files
set(SRCS
  Hit.cxx
)

set(HEADERS
  Hit.h
)

Set(LINKDEF DataFormatsLinkDef.h)
Set(LIBRARY_NAME DataFormats)
Set(DEPENDENCIES)
Set(DictName "G__${LinName}Dict.cxx")
generate_dictionary("DataFormats" "${DictName}" "${HEADERS}" "${INCLUDE_DIRECTORIES}")
generate_rootmap("DataFormats" "${INCLUDE_DIRECTORIES}" "${CMAKE_CURRENT_SOURCE_DIR}/${DataFormats}LinkDef.h")
GENERATE_LIBRARY()
