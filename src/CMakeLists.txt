cmake_minimum_required(VERSION 2.8)

project(nbody)

add_executable(brutforce-serial brutforce_serial.cpp)
add_executable(brutforce-parallel brutforce_parallel.cpp)
add_executable(quadtree-serial quadtree_serial.cpp Quadtree.cpp)
add_executable(quadtree-parallel quadtree_parallel.cpp Quadtree.cpp)

# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILER /usr/bin/mpicxx)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
set(CMAKE_CXX_FLAGS:STRING="-std=c++11")
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(brutforce-serial ${MPI_LIBRARIES})
target_link_libraries(brutforce-parallel ${MPI_LIBRARIES})
target_link_libraries(quadtree-serial ${MPI_LIBRARIES})
target_link_libraries(quadtree-parallel ${MPI_LIBRARIES})