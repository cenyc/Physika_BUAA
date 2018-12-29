cmake_minimum_required(VERSION 3.10)
project(Physika_CloudRender)
include("$ENV{CMAKE_CUDA}")

set(CMAKE_CXX_STANDARD 11)
#include_directories(Physika_Arrays)
#include_directories(Physika_Utilities)
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}")
#set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -m64")
#set(CUDA_64_BIT_DEVICE_CODE ON)
#set(BUILD_USE_64BITS on)
#set(NVC)
#SET(CMAKE_BUILD_TYPE Release)

#czn
message(STATUS  "=============cmake Physika_CloudRender========")

message(STATUS "cmake cuda path: " "$ENV{CMAKE_CUDA}" )
include_directories(.)#;../Cloud_Dependency/trimesh2/include;../Cloud_Dependency/SDFGen
aux_source_directory(. srcs)
#message(STATUS  "=============" ${srcs})

add_executable(Physika_CloudRender ${srcs})
add_library(Physika_CloudRender_lib ${srcs})
target_link_libraries(Physika_CloudRender Physika_CloudDepdency_lib SDFGen_lib)

