#cmake_minimum_required (VERSION 2.6)
#project (TSPEED)

FILE(GLOB components "src/*.cpp")

add_library(tspeed ${components})

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX
	"/usr/local/" CACHE PATH "/usr/local install prefix" FORCE
	)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
install(TARGETS tspeed DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
FILE(GLOB headers "include/*.hpp")
install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/tspeed/)

message(STATUS "Installing to ${CMAKE_INSTALL_PREFIX}, to change rerun cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/dir")

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Default build type: Release. Use -DCMAKE_BUILD_TYPE=Debug build with debug symbols.")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")
find_package(Eigen3 REQUIRED)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_LIBCXX11)
if(NOT HAS_LIBCXX11)
    message(FATAL_ERROR "Current compiler does not suppport -std=c++11 option. Exiting Cmake.")
endif()




include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/include", "${EIGEN3_INCLUDE_DIR}")

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -msse4 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall")


