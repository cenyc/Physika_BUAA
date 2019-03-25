cmake_minimum_required(VERSION 3.10)
project(Cloud_Gen_Based_Euler)

set(CMAKE_CXX_STANDARD 11)

# add_executable(Cloud_Gen_Based_Euler demo_test.cpp Cloud_gen.cpp)

aux_source_directory(. Cloud_Gen_Based_Euler)
add_executable(Cloud_Gen_Based_Euler ${Cloud_Gen_Based_Euler})

# find_package(OpenGL REQUIRED)
# find_package(GLUT REQUIRED)

# include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})
# target_link_libraries(Cloud_Gen_Based_Euler ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})

add_library(Cloud_Gen_Based_Euler_lib ${Cloud_Gen_Based_Euler})
