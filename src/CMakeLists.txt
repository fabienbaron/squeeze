cmake_minimum_required(VERSION 3.0)
project(squeeze C)

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

set(CMAKE_BUILD_TYPE release)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -march=native -mfpmath=sse -ftree-vectorize -ggdb -fno-omit-frame-pointer")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fno-omit-frame-pointer -O2") #Note: valgrind chokes on march=native for some platforms

file(GLOB SOURCE *.h squeeze.c free_fits.c read_fits.c chkmalloc.c ../lib/oifitslib/src/oifitslib/exchange.h ../lib/oifitslib/src/oifitslib/free_fits.c ../lib/oifitslib/src/oifitslib/read_fits.c ../lib/oifitslib/src/oifitslib/chkmalloc.h ../lib/oifitslib/src/oifitslib/chkmalloc.c ../lib/oifitslib/src/oifitslib/alloc_fits.c )

# Now add the binary
add_executable(squeeze ${SOURCE})

target_link_libraries(squeeze m pthread ${CMAKE_CURRENT_SOURCE_DIR}/../lib/cfitsio/libcfitsio.a ${CMAKE_CURRENT_SOURCE_DIR}/../lib/rngstreams/src/.libs/librngstreams.a)   


