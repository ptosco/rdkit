# - Try to find OpenMM lib
# Once done this will define
#
#  OPENMM_FOUND - system has OpenMM lib
#  OPENMM_INCLUDE_DIR - the OpenMM include directory
#  OPENMM_LIBRARIES - the OpenMM libraries we want to link

if(OPENMM_INCLUDE_DIR AND OPENMM_LIBRARIES)
  # in cache already
  set(OPENMM_FOUND TRUE)
else()
  find_path(OPENMM_INCLUDE_DIR NAMES OpenMM.h PATHS /usr/local/openmm/include)
  foreach(openMMLibraryName OpenMM OpenMMMMFF)
    find_library(OPENMM_LIBRARY NAMES ${openMMLibraryName})
    if((${OPENMM_LIBRARY} STREQUAL "OPENMM_LIBRARY-NOTFOUND") AND OPENMM_LIBRARY_DIR)
      find_library(OPENMM_LIBRARY PATHS ${OPENMM_LIBRARY_DIR} NAMES ${openMMLibraryName})
    endif()
    if(NOT (${OPENMM_LIBRARY} STREQUAL "OPENMM_LIBRARY-NOTFOUND"))
      if(OPENMM_LIBRARIES)
        set(OPENMM_LIBRARIES "${OPENMM_LIBRARIES};")
      endif()
      set(OPENMM_LIBRARIES "${OPENMM_LIBRARIES}${OPENMM_LIBRARY}")
    endif()
    unset(OPENMM_LIBRARY CACHE)
  endforeach()
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(OpenMM
	  "Could NOT find OpenMM in system locations" OPENMM_LIBRARIES
	  OPENMM_INCLUDE_DIR)
  mark_as_advanced(OPENMM_INCLUDE_DIR OPENMM_LIBRARIES)
  message("OPENMM_LIBRARIES = ${OPENMM_LIBRARIES}")
endif()
