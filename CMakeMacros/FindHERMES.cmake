INCLUDE(FindPackageHandleStandardArgs)

if(${HERMES_VERSION} STREQUAL "HERMES2D")

IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
  FIND_LIBRARY(HERMES_LIBRARY NAMES libhermes2d-debug hermes2d-debug PATHS ${HERMES_DIRECTORY}/hermes2d /usr/lib /usr/local/lib) #Only debug!
  MESSAGE("debug mode")
ELSE(CMAKE_BUILD_TYPE MATCHES DEBUG)
 FIND_LIBRARY(HERMES_LIBRARY NAMES libhermes2d hermes2d libhermes2d-debug hermes2d-debug PATHS ${HERMES_DIRECTORY}/hermes2d /usr/lib /usr/local/lib)
ENDIF() 

  MESSAGE(${HERMES_LIBRARY})
  FIND_PATH(HERMES_INCLUDE hermes2d.h ${HERMES_DIRECTORY}/hermes2d/include /usr/include/hermes2d /usr/local/include/hermes2d)
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(HERMES DEFAULT_MSG HERMES_LIBRARY HERMES_INCLUDE)
endif(${HERMES_VERSION} STREQUAL "HERMES2D")
