cmake_minimum_required(VERSION 2.8) 
project(squeeze)

# If you want to override compiler please uncomment the line below
#SET(CMAKE_C_COMPILER pathto/mycompiler)

# Compile with C11 support.
include(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG("-std=c11" COMPILER_SUPPORTS_C11)
if(COMPILER_SUPPORTS_C11)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
else()
    message(STATUS "The compiler ${CMAKE_C_COMPILER} has no C11 support. Please use a different C compiler.")
endif()

# Compile with OpenMP support
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
else()
    message(STATUS "The compiler ${CMAKE_C_COMPILER} has no OpenMP support. Please use a different C compiler if you want to use OpenMP capabilities.")
endif()

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS_DEBUG "-Ofast -march=native -mfpmath=sse -flto -ftree-vectorize -ggdb -fno-omit-frame-pointer")

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
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/cfitsio
  CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/lib/cfitsio/configure --silent --enable-reentrant --enable-ssse3
#  PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/lib/cfitsio
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND ""
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/cfitsio)
set(LIBS ${LIBS} "${CMAKE_CURRENT_SOURCE_DIR}/lib/cfitsio/libcfitsio.a")

# Build the main directory, always
add_subdirectory(src)


