project(vegas C)

cmake_minimum_required(VERSION 2.8)

find_package(MPI)

include_directories(${MPI_C_INCLUDE_PATH})
set(CMAKE_COMPILE_FLAGS ${CMAKE_COMPILE_FLAGS} ${MPI_C_COMPILE_FLAGS})

add_executable(vegas_test main.c vegas.c vegas.h matrix.h)
target_link_libraries(vegas_test ${MPI_C_LIBRARIES})