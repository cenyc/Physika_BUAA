cmake_minimum_required(VERSION 3.10)
project(Physika_CloudDepdency)

set(CMAKE_CXX_STANDARD 11)

include_directories(./trimesh2/include)
set(srcs "")
aux_source_directory(./trimesh2/libsrc srcs)
aux_source_directory(./trimesh2/utilsrc srcs)
add_library(Physika_CloudDepdency_lib  ${srcs})
add_subdirectory(./SDFGen)
#add_executable(Physika_CloudRende88r ./renderMain.cpp )
#target_link_libraries(Physika_CloudDepdency_lib trimesh2_lib)
