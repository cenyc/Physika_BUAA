cmake_minimum_required(VERSION 3.10)
project(Physika)

set(CMAKE_CXX_STANDARD 11)
include_directories(Physika_Src)


#message(STATUS "envs:")
#message(STATUS "CUDA_API_PATH $ENV{CUDA_API_PATH}")


set(CMAKE_BUILD_TYPE debug)


#set(PROJECT_BINARY_DIR   ${CMAKE_SOURCE_DIR}/cmake_build)
set(CMAKE_BINARY_DIR    ${CMAKE_SOURCE_DIR}/cmake_build)
set(CMAKE_CURRENT_BINARY_DIR     ${CMAKE_SOURCE_DIR}/cmake_build)

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#message(STATUS "++++++++++++" ${CMAKE_BINARY_DIR} " " ${PROJECT_BINARY_DIR} " " ${CMAKE_SOURCE_DIR} ${CMAKE_HOME_DIRECTORY})

add_subdirectory(Physika_Src)

add_executable(Physika main.cpp )
#file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
