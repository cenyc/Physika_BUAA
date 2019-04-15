cmake_minimum_required(VERSION 3.10)
project(Cloud_HalfMesh)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS -m64)
# add_executable(Cloud_HalfMesh main.cpp Cloud.cpp CloudInterior.cpp image.cpp Mesh.cpp Perlin.cpp Pixel.cpp sky.cpp Sun.cpp Tool.cpp)

aux_source_directory(. Cloud_HalfMesh)

set(OpenCV_DIR D:/cv3.0.0/opencv/build)
find_package(OpenCV REQUIRED)

# include_directories(${OpenCV_INCLUDE_DIRS})
include_directories("D:/cv3.0.0/opencv/build/include")

# add_library(Cloud_HalfMesh_lib ${Cloud_HalfMesh})

add_executable(Cloud_HalfMesh ${Cloud_HalfMesh})

# target_link_directories(Cloud_HalfMesh ${OpenCV_LIBS})