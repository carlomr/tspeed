#cmake_minimum_required (VERSION 2.6)
#project (TSPEED)

link_directories(${CMAKE_INSTALL_PREFIX}/lib)
add_executable(Wedge src/wedge.cpp)
add_executable(Lamb src/Lamb.cpp)

install(TARGETS Wedge Lamb DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}  )
target_link_libraries(Wedge tspeed)
target_link_libraries(Lamb tspeed)



include_directories ("${CMAKE_SOURCE_DIR}/lib/include", ${EIGEN3_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -msse2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall")



