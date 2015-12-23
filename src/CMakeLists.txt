cmake_minimum_required(VERSION 2.8)

project(nbody)

find_package(MPI REQUIRED)

set(CMAKE_CXX_FLAGS "-std=c++0x" CACHE STRING "Compilation flags")
include_directories(${MPI_CXX_INCLUDE_PATH})

set(_targets brutforce-serial brutforce-parallel quadtree-serial quadtree-parallel)

add_executable(brutforce-serial brutforce_serial.cpp)
add_executable(brutforce-parallel brutforce_parallel.cpp)
add_executable(quadtree-serial quadtree_serial.cpp Quadtree.cpp)
add_executable(quadtree-parallel quadtree_parallel.cpp Quadtree.cpp)

foreach(_target ${_targets})
  target_link_libraries(${_target} ${MPI_CXX_LIBRARIES})
endforeach()
