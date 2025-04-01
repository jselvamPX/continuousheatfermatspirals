cmake_minimum_required(VERSION 3.16)
project(Fermat3D)
set(CMAKE_CXX_STANDARD 17)

# Include and fetch Eigen3
include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)

# Fetch and configure OpenMesh
FetchContent_Declare(
    openmesh
    GIT_REPOSITORY https://www.graphics.rwth-aachen.de:9000/OpenMesh/OpenMesh.git
    GIT_TAG OpenMesh-11.0
    GIT_SHALLOW TRUE
)

# Configure OpenMesh build options before making it available
set(BUILD_APPS OFF CACHE BOOL "Build OpenMesh apps")
set(OPENMESH_DOCS OFF CACHE BOOL "Build OpenMesh documentation")
set(BUILD_TESTING OFF CACHE BOOL "Build OpenMesh tests")

FetchContent_MakeAvailable(openmesh)

# Get OpenMesh source and binary directories
FetchContent_GetProperties(openmesh)

# Fetch libigl
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(libigl)

# Add your executable
add_executable(heatisoline 
    HeatFermat3D.cpp 
    Spirals.cpp 
    SpiralDt.cpp
)

# Include directories for OpenMesh
target_include_directories(heatisoline PRIVATE 
    ${openmesh_SOURCE_DIR}/src
    ${openmesh_BINARY_DIR}/src
)

# Link dependencies
target_link_libraries(heatisoline PRIVATE
    igl::core
    OpenMeshCore
    OpenMeshTools
    Eigen3::Eigen
)