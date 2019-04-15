PROJECT(MeshDeform)

cmake_minimum_required(VERSION 3.12)

#link_libraries("${CMAKE_SOURCE_DIR}/build/Physika_Src/libPhysika_Core.a")

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

find_package(Eigen3 REQUIRED)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
find_package(OpenMesh REQUIRED)

include_directories(
	${OPENMESH_INCLUDE_DIR}
	${EIGEN3_INCLUDE_DIR}
	${PROJECT_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/Physika_Src
	)

#link_directories(
#	${OPENMESH_LIBRARY_DIR}
#	)

target_link_libraries(${PROJECT_NAME} ${OPENMESH_LIBRARY_DIR}/libOpenMeshCore.so ${OPENMESH_LIBRARY_DIR}/libOpenMeshTools.so)
