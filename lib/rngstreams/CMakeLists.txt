cmake_minimum_required(VERSION 3.0) 
project(librngstream)

set(CMAKE_BUILD_TYPE Release)
set(CONFIGURE_COMMAND ./configure )

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

# Include and build Google Test framework:
#add_subdirectory(src)
