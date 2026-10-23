cmake_minimum_required(VERSION 3.10)
project(Simulation)

find_package(ROOT CONFIG REQUIRED COMPONENTS Core EG Geom RIO Tree)

find_package(Geant4 REQUIRED)

find_package(VMC REQUIRED)

find_package(VGM REQUIRED)

find_package(Geant4VMC REQUIRED)

#Create a library called "libSimulation"
Set(LibName "Simulation")

get_target_property(ROOT_INCLUDE_DIRS ROOT::Core INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(ROOT_LIBRARY_DIR ROOT::Core INTERFACE_LINK_DIRECTORIES)

set(INCLUDE_DIRECTORIES 
    ${ROOT_INCLUDE_DIRS}
    ${PYTHIA_INCLUDE_DIR}/.. 
    ${Geant4_INCLUDE_DIRS} 
    ${Geant4VMC_INCLUDE_DIRS}
    ${VMC_INCLUDE_DIRS}
    ${VGM_INCLUDE_DIRS}
    ${PHOS256_SOURCE_DIR}/utils
    ${PHOS256_SOURCE_DIR}/dataformats
    ${PHOS256_SOURCE_DIR}/raw
    ${CMAKE_CURRENT_LIST_DIR}
)

set(LINK_DIRECTORIES 
    ${ROOT_LIBRARY_DIR}                       
    ${PYTHIA_LIBRARY_SHARED}
    ${Geant4VMC_LIBRARY_DIR}
)

link_directories(${LINK_DIRECTORIES})

#List of source files
set(SRCS  Hall.cxx Magnet.cxx Phos.cxx Stack.cxx GenBox.cxx GenPythia.cxx GenPolPhoton.cxx GenMuon.cxx MagField.cxx Digitizer.cxx Simulation.cxx RunConfigurationPol.cxx)

set(HDRS Hall.h Magnet.h Phos.h Stack.h GenBox.h GenPythia.h GenPolPhoton.h GenMuon.h MagField.h Digitizer.h Simulation.h RunConfigurationPol.h)

#Generate the dictionary
get_directory_property(incdirs INCLUDE_DIRECTORIES)
generate_dictionary("${LibName}" "${LibName}LinkDef.h" "${HDRS}" "${INCLUDE_DIRECTORIES}")

set(ROOT_LIBRARIES_NAMES Core EG Geom RIO Tree VMC)
set(LIBDEPS ${ROOT_LIBRARIES_NAMES} Geant4 Geant4VMC VGM)
generate_rootmap("${LibName}" "${LIBDEPS}" "${CMAKE_CURRENT_SOURCE_DIR}/${LibName}LinkDef.h")

#Create an object to be reused
add_library(${LibName}-object OBJECT ${SRCS} G__${LibName}.cxx)

#Add a library to the project using the object
add_library(${LibName} SHARED $<TARGET_OBJECTS:${LibName}-object>)

target_link_libraries(${LibName} 
    ROOT::Core
    ROOT::EG
    ROOT::Geom
    ROOT::RIO
    ROOT::Tree
    ${Geant4_LIBRARIES}
    ${Geant4VMC_LIBRARIES}
    ${VGM_LIBRARIES}
)


#Setting the correct headers for the object
target_include_directories(${LibName}-object PUBLIC $<TARGET_PROPERTY:${LibName},INCLUDE_DIRECTORIES>)
set_target_properties(${LibName}-object PROPERTIES COMPILE_DEFINITIONS $<TARGET_PROPERTY:${LibName},COMPILE_DEFINITIONS>)

#Public include folders
target_include_directories(${LibName} PUBLIC ${INCLUDE_DIRECTORIES})

set(LibName_COMPILE_FLAGS)
set(LibName_LINK_FLAGS)

set_target_properties(${LibName}-object PROPERTIES COMPILE_FLAGS "${LibName_COMPILE_FLAGS}")
set_target_properties(${LibName} PROPERTIES LINK_FLAGS "${LibName_LINK_FLAGS}")

#Installation
install(TARGETS ${LibName} ARCHIVE DESTINATION lib LIBRARY DESTINATION lib)
install(FILES ${HDRS} DESTINATION include)