# - Try to find OpenMM lib
# Once done this will define
#
#  OPENMM_FOUND - system has OpenMM lib
#  OPENMM_INCLUDE_DIRS - the OpenMM include directory
#  OPENMM_LIBRARIES - the OpenMM libraries we want to link

if(OPENMM_INCLUDE_DIR AND OPENMM_LIBRARY)
  # in cache already
  set(OPENMM_FOUND TRUE)
else()
  find_path(OPENMM_INCLUDE_DIR NAMES OpenMM.h PATHS /usr/local/openmm/include)
  find_library(OPENMM_LIBRARY NAMES OpenMM)
  if ((${OPENMM_LIBRARY} STREQUAL "OPENMM_LIBRARY-NOTFOUND") AND OPENMM_LIBRARY_DIR)
    find_library(OPENMM_LIBRARY PATHS ${OPENMM_LIBRARY_DIR} NAMES OpenMM)
  endif()
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(OPENMM
	  "Could NOT find OpenMM in system locations" OPENMM_LIBRARY
	  OPENMM_INCLUDE_DIR)
  set(OPENMM_INCLUDE_DIRS ${OPENMM_INCLUDE_DIR})
  set(OPENMM_LIBRARIES ${OPENMM_LIBRARY})
  mark_as_advanced(OPENMM_INCLUDE_DIR OPENMM_LIBRARY)
endif()
