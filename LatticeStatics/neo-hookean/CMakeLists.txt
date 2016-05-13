CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.3 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../../ ../../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

OPTION(CREATE_LIB "Create a library, otherwise create executable" OFF)

PROJECT(neo-hookean)

IF(CREATE_LIB)
  ADD_LIBRARY(neo-hookean SHARED # or STATIC
    neo-hookean.cpp
    )
  ADD_DEFINITIONS(-DCREATE_LIBRARY)
  DEAL_II_SETUP_TARGET(neo-hookean)
ELSE(CREATE_LIB)
  SET(TARGET "neo-hookean")
  SET(TARGET_SRC ${TARGET}.cpp)
  PROJECT(${TARGET})
  DEAL_II_INVOKE_AUTOPILOT()
ENDIF(CREATE_LIB)
