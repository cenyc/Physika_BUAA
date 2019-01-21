cmake_minimum_required(VERSION 3.10)
project(trimesh2)

set(CMAKE_CXX_STANDARD 11)

include_directories(./include)
set(srcs "")
aux_source_directory(./libsrc srcs)
aux_source_directory(./utilsrc srcs)

add_library(trimesh2_lib   ${srcs})
#add_executable(Physika_CloudRende88r ./renderMain.cpp )
#target_link_libraries(Physika_Dependency Physika_Dependency_lib)
