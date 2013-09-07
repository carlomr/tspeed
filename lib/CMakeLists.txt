#cmake_minimum_required (VERSION 2.6)
#project (TSPEED)

FILE(GLOB components "src/*.cpp")

add_library(tspeed ${components})

install(TARGETS tspeed DESTINATION ${PROJECT_SOURCE_DIR})

message("${PROJECT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")
find_package(Eigen3 REQUIRED)
message("${CMAKE_MODULE_PATH}")

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_LIBCXX11)
if(NOT HAS_LIBCXX11)
    message(FATAL_ERROR "Current compiler does not suppport -std=c++11 option. Exiting Cmake.")
endif()




include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/include", "${EIGEN3_INCLUDE_DIR}")

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -msse4 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall")


