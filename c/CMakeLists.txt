cmake_minimum_required(VERSION 3.10)
project(Smart_Antennas)
set(CMAKE_C_STANDARD 11)

find_package(osqp REQUIRED)

include_directories("include")
include_directories("src")

set(SOURCE_FILES main.c include/core.h src/core.c)

add_executable(Smart_Antennas ${SOURCE_FILES})

target_link_libraries(Smart_Antennas PRIVATE osqp::osqp)
