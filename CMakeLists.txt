cmake_minimum_required (VERSION 2.6)
project (TSPEED)
#add_executable(Wedge wedge.cpp)

#install(TARGETS Wedge DESTINATION ${PROJECT_SOURCE_DIR}  )
#target_link_libraries(Wedge tspeed)

add_subdirectory(lib)
add_subdirectory(Examples)
add_subdirectory(doc)

#install(TARGETS tspeed DESTINATION ${PROJECT_SOURCE_DIR})


set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(Eigen3 REQUIRED)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" HAS_LIBCXX11)
if(NOT HAS_LIBCXX11)
    message(FATAL_ERROR "Current compiler does not suppport -std=c++11 option. Exiting Cmake.")
endif()

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Default build type: Release. Use -DCMAKE_BUILD_TYPE=Debug build with debug symbols.")
    set(CMAKE_BUILD_TYPE "Release")
endif()

message(STATUS "Run make doc to build Doxygen documentation")
include_directories ("${PROJECT_SOURCE_DIR}/lib/include", ${EIGEN3_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -msse4 ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall ")



