cmake_minimum_required(VERSION 3.5)
project(Fermat3D)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(LIBIGL REQUIRED QUIET)

find_package(Eigen3 CONFIG REQUIRED)

find_package(OpenMesh CONFIG REQUIRED)

add_executable(heatisoline HeatFermat3D.cpp Spirals.cpp SpiralDt.cpp)

target_include_directories(heatisoline PRIVATE ${OPENMESH_INCLUDE_DIR} ${LIBIGL_INCLUDE_DIR})

target_link_libraries(heatisoline PRIVATE
        Eigen3::Eigen
        igl
        ${OPENMESH_LIBRARIES} )