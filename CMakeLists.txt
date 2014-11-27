cmake_minimum_required(VERSION 2.8) 
project(squeeze)

# If you want to override compiler please uncomment the line below
#SET(CMAKE_C_COMPILER pathto/mycompiler)

set(CMAKE_BUILD_TYPE release)

# Set some CMake properties:
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)
#include_directories(${CMAKE_SOURCE_DIR}/include)

# We need to compile the following libraries

# rngstreams
include(ExternalProject)
ExternalProject_Add(
  librngstreams
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/rngstreams
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/lib/rngstreams/configure --silent
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND ""
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/rngstreams)
set(LIBS ${LIBS} "${CMAKE_CURRENT_SOURCE_DIR}/lib/rngstreams/librngstreams.a")

# cfitsio
ExternalProject_Add(
  libcfitsio
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/cfitsio
  CMAKE_ARGS "-DUSE_PTHREADS=ON"
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
#  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/lib/cfitsio/configure --silent --enable-reentrant --enable-ssse3
#  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib/cfitsio 
  INSTALL_COMMAND ""
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/cfitsio)
set(LIBS ${LIBS} "${CMAKE_CURRENT_SOURCE_DIR}/lib/cfitsio/libcfitsio.a")

# Build the main directory, always
add_subdirectory(src)


