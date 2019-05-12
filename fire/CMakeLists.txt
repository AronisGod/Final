cmake_minimum_required(VERSION 3.14)
project(fire)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(fire
        Fire.cpp
        Fire.h
        main.cpp)
