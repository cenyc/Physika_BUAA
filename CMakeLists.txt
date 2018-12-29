cmake_minimum_required(VERSION 3.10)
project(Physika)

set(CMAKE_CXX_STANDARD 11)
include_directories(Physika_Src)



#czn


message(STATUS "envs:")
message(STATUS "CUDA_API_PATH $ENV{CUDA_API_PATH}")


add_subdirectory(Physika_Src/Physika_Cloud/Cloud_Dependency)

ADD_SUBDIRECTORY(Physika_Src/Physika_Cloud/Cloud_Render)
add_executable(Physika main.cpp )
#target_link_libraries(Physika Physika_Core_lib)
#target_link_libraries(Physika Physika_Dependency_lib)
