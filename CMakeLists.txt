CMAKE_MINIMUM_REQUIRED( VERSION 2.8.7 )

SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -O3 -fopenmp")
PROJECT(fluid_simulation)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${fluid_simulation_SOURCE_DIR}/bin)

SET(PROJECT_HEADER
        ${PROJECT_SOURCE_DIR}/include/CIsoSurface.h
        ${PROJECT_SOURCE_DIR}/include/Vectors.h
        ${PROJECT_SOURCE_DIR}/include/stdafx.h
        ${PROJECT_SOURCE_DIR}/include/display.h
        ${PROJECT_SOURCE_DIR}/include/fluidCube.h
        ${PROJECT_SOURCE_DIR}/include/fluidCube2D.h
        ${PROJECT_SOURCE_DIR}/include/fluidCube3D.h
        ${PROJECT_SOURCE_DIR}/include/Configer.h
        ${PROJECT_SOURCE_DIR}/include/Logger.h
        ${PROJECT_SOURCE_DIR}/include/trianglemesh.h
        ${PROJECT_SOURCE_DIR}/include/quadmesh.h)

SET(PROJECT_SRC
        ${PROJECT_SOURCE_DIR}/src/CIsoSurface.cpp
        ${PROJECT_SOURCE_DIR}/src/Vectors.cpp
        ${PROJECT_SOURCE_DIR}/src/display.cpp
        ${PROJECT_SOURCE_DIR}/src/fluidCube.cpp
        ${PROJECT_SOURCE_DIR}/src/fluidCube2D.cpp
        ${PROJECT_SOURCE_DIR}/src/fluidCube3D.cpp
        ${PROJECT_SOURCE_DIR}/src/Configer.cpp
        ${PROJECT_SOURCE_DIR}/src/Logger.cpp
        ${PROJECT_SOURCE_DIR}/src/quadmesh.cpp
        ${PROJECT_SOURCE_DIR}/src/trianglemesh.cpp
        ${PROJECT_SOURCE_DIR}/src/main.cpp)

INCLUDE_DIRECTORIES(include)

ADD_EXECUTABLE(fluid_simulation ${PROJECT_HEADER} ${PROJECT_SRC})

# Eigen
FIND_PACKAGE(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (fluid_simulation Eigen3::Eigen)

# opencv
FIND_PACKAGE(OpenCV REQUIRED)
TARGET_LINK_LIBRARIES(fluid_simulation ${OpenCV_LIBS})

# opengl
FIND_PACKAGE(OpenGL REQUIRED)
LIST(APPEND LIBS ${OPENGL_gl_LIBRARY})
INCLUDE_DIRECTORIES(${OPENGL_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(fluid_simulation ${OPENGL_LIBRARIES})

# GLUT
FIND_PACKAGE(GLUT REQUIRED)
INCLUDE_DIRECTORIES(${GLUT_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(fluid_simulation ${GLUT_LIBRARIES})
