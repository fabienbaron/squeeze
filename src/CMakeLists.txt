cmake_minimum_required(VERSION 2.8)
project(squeeze C)

file(GLOB SOURCE *.h squeeze.c free_fits.c read_fits.c)

# Now add the binary
add_executable(squeeze ${SOURCE})

target_link_libraries(squeeze m pthread ${CMAKE_CURRENT_SOURCE_DIR}/../lib/cfitsio/libcfitsio.a ${CMAKE_CURRENT_SOURCE_DIR}/../lib/rngstreams/src/.libs/librngstreams.a)   

