# Create a library called "libDataFormats" which includes the source files given in the array .
Set(LibName "DataFormats")

set(INCLUDE_DIRECTORIES
  ${ROOT_INCLUDE_DIRS}
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
set(SRCS Hit.cxx Digit.cxx
)

set(HDRS Hit.h Digit.h
)

# Generate the dictionary
# It will create G_ARG1.cxx and G_ARG1.h / ARG1 = function first argument
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${LibName}" "${LibName}LinkDef.h" "${HDRS}" "${INCLUDE_DIRECTORIES}")

set(ROOT_DEPENDENCIES Core)

# Generate the ROOT map
# Dependecies
set(LIBDEPS ${ROOT_DEPENDENCIES})
generate_rootmap("${LibName}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${LibName}LinkDef.h")

# Create an object to be reused in case of static libraries
# Otherwise the sources will be compiled twice
add_library(${LibName}-object OBJECT ${SRCS} G__${LibName}.cxx)
# Add a library to the project using the object
add_library(${LibName} SHARED $<TARGET_OBJECTS:${LibName}-object>)
target_link_libraries(${LibName} ${ROOT_DEPENDENCIES})

# Setting the correct headers for the object as gathered from the dependencies
target_include_directories(${LibName}-object PUBLIC $<TARGET_PROPERTY:${LibName},INCLUDE_DIRECTORIES>)
set_target_properties(${LibName}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${LibName},COMPILE_DEFINITIONS>)

# Public include folders that will be propagated to the dependecies
target_include_directories(${LibName} PUBLIC ${INCLUDE_DIRECTORIES})

set(LibName_COMPILE_FLAGS)
set(LibName_LINK_FLAGS)

# Setting compilation flags for the object
set_target_properties(${LibName}-object PROPERTIES COMPILE_FLAGS "${LibName_COMPILE_FLAGS}")
# Setting the linking flags for the library
set_target_properties(${LibName} PROPERTIES LINK_FLAGS "${LibName_LINK_FLAGS}")

# Installation
install(TARGETS ${LibName}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES ${HDRS} DESTINATION include)
