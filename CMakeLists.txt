project(ATMesh)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
cmake_minimum_required(VERSION 2.6)

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})
include_directories(${OPENGL_INCLUDE_DIR}/Headers)
link_directories(${OPENGL_LIBRARY_DIRS})

find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIR} ${GLUT_INCLUDE_DIR}/GL ${GLUT_INCLUDE_DIR}/GLUT)
link_directories(${GLUT_LIBRARY_DIRS})

find_package(GLEW REQUIRED)
include_directories(${GLEW_INCLUDE_DIRS})
include_directories(${GLEW_INCLUDE_DIRS}/GL)
link_directories(${GLEW_LIBRARY_DIRS})

find_package(SuiteSparse REQUIRED)
include_directories(${SUITESPARSE_INCLUDE_DIRS})

find_package(Umfpack REQUIRED)
include_directories(${UMFPACK_INCLUDES})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


find_package(Boost 1.46.0 REQUIRED COMPONENTS program_options)
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories(${Boost_LIBRARY_DIRS})

add_compile_options("-std=c++11")

add_library(dec STATIC
    DenseMatrix.cpp
    Shader.cpp
    Variable.cpp
    Real.cpp
    LinearContext.cpp
    Edge.cpp
    LinearPolynomial.cpp
    MeshIO.cpp
    Vertex.cpp
    Face.cpp
    Mesh.cpp
    LinearEquation.cpp
    Complex.cpp
    HalfEdge.cpp
    Camera.cpp
    LinearSystem.cpp
    SparseMatrix.cpp
    Vector.cpp
    Image.cpp
    Quaternion.cpp
    )

target_link_libraries(dec
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    ${GLEW_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${UMFPACK_LIBRARIES}
    ${Boost_LIBRARIES}
    )

include_directories(".")
add_subdirectory(positions)
