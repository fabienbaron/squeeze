cmake_minimum_required(VERSION 3.20)

project(squeeze
  VERSION 2.8.0
  DESCRIPTION "Image reconstruction for optical interferometry by MCMC"
  HOMEPAGE_URL "https://github.com/fabienbaron/squeeze"
  LANGUAGES C)

# ---------------------------------------------------------------------------
# Language / build configuration
# ---------------------------------------------------------------------------
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release RelWithDebInfo MinSizeRel)
endif()

# Defaults to ./bin (as the docs and helper scripts expect) but stays
# overridable with -DCMAKE_RUNTIME_OUTPUT_DIRECTORY=... so that separate build
# trees do not clobber each other's binary.
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
endif()

option(SQUEEZE_NATIVE_ARCH
  "Optimise for this machine's CPU (-march=native). Turn OFF for binaries that must run on other hardware."
  ON)

# In-source builds leave CMakeCache.txt/CMakeFiles littering the repo root.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR
    "In-source builds are not supported. Run:\n"
    "  cmake -B build && cmake --build build -j")
endif()

# ---------------------------------------------------------------------------
# OpenMP (parallel tempering / simulated annealing)
# ---------------------------------------------------------------------------
find_package(OpenMP COMPONENTS C)
if(NOT OpenMP_C_FOUND)
  message(WARNING
    "No OpenMP support found for ${CMAKE_C_COMPILER_ID}. SQUEEZE will build, but "
    "parallel tempering and parallel annealing will run on a single thread.")
endif()

# ---------------------------------------------------------------------------
# CFITSIO (vendored: lib/cfitsio)
# ---------------------------------------------------------------------------
# Force-set as cache entries so they win over CFITSIO's own option() defaults.
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)  # static: avoids runtime lookup issues
set(USE_CURL          OFF CACHE BOOL "" FORCE)  # no remote FITS access needed
set(USE_PTHREADS      ON  CACHE BOOL "" FORCE)  # SQUEEZE reads FITS from OpenMP threads
set(TESTS             OFF CACHE BOOL "" FORCE)  # skip TestProg/cookbook
set(UTILS             OFF CACHE BOOL "" FORCE)  # skip fpack/funpack/fitsverify/...
set(ITERPROGS         OFF CACHE BOOL "" FORCE)

add_subdirectory(lib/cfitsio EXCLUDE_FROM_ALL)

# CFITSIO's CMakeLists never exports its own headers to consumers, so `#include
# <fitsio.h>` would not resolve through the target alone. BUILD_INTERFACE keeps
# the path out of CFITSIO's own install(EXPORT) rules, which would otherwise
# reject a source-tree path.
target_include_directories(cfitsio
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/lib/cfitsio>)

# ---------------------------------------------------------------------------
# RngStreams (vendored: lib/rngstreams)
# ---------------------------------------------------------------------------
# Upstream ships an autotools build; the library is two files, so compile it
# directly rather than driving ./configure through ExternalProject.
add_library(rngstreams STATIC lib/rngstreams/src/RngStream.c)
target_include_directories(rngstreams PUBLIC "${PROJECT_SOURCE_DIR}/lib/rngstreams/src")
target_link_libraries(rngstreams PRIVATE m)

# ---------------------------------------------------------------------------
# OIFITSlib exchange layer (submodule: lib/oifitslib)
# ---------------------------------------------------------------------------
set(OITABLE_DIR "${PROJECT_SOURCE_DIR}/lib/oifitslib/src/oifitslib")
if(NOT EXISTS "${OITABLE_DIR}/exchange.h")
  message(FATAL_ERROR
    "The oifitslib submodule is missing. Run:\n"
    "  git submodule update --init")
endif()

# Mirrors upstream's `oitable` target: the low-level table layer only. The
# higher-level `oifits` target needs GLib, which SQUEEZE does not use.
add_library(oitable STATIC
  "${OITABLE_DIR}/read_fits.c"
  "${OITABLE_DIR}/write_fits.c"
  "${OITABLE_DIR}/alloc_fits.c"
  "${OITABLE_DIR}/free_fits.c"
  "${OITABLE_DIR}/chkmalloc.c")
target_include_directories(oitable PUBLIC "${OITABLE_DIR}")
target_link_libraries(oitable PUBLIC cfitsio PRIVATE m)

# ---------------------------------------------------------------------------
add_subdirectory(src)
