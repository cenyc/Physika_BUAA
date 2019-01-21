cmake_minimum_required(VERSION 3.12)
project(Render_Harris)

set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES main.cpp)

include_directories(.)

link_directories(${LINK_DIR})

#add_executable(Render_Harris
#        Camera.cpp
#        Camera.h
#        Color.h
#        Common.h
#        Frustum.cpp
#        Frustum.h
#        GLText.cpp
#        GLText.h
#        main.cpp
#        Matrix.cpp
#        Matrix.h
#        SkyDome.cpp
#        SkyDome.h
#        Vector.h
#        VolumeClouds.cpp
#        VolumeClouds.h
#        utils.h
#        )

aux_source_directory(. temp)
add_executable(Render_Harris ${temp})
include_directories(/usr/include/IL)
target_link_libraries(${PROJECT_NAME} /usr/lib/x86_64-linux-gnu/libIL.so /usr/lib/x86_64-linux-gnu/libILU.so /usr/lib/x86_64-linux-gnu/libILUT.so)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(SDL REQUIRED)

include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${SDL_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${SDL_LIBRARY})
