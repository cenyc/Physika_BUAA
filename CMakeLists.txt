cmake_minimum_required(VERSION 3.10)
project(Physika)

set(CMAKE_CXX_STANDARD 11)
include_directories(Physika_Src)
add_executable(Physika main.cpp)