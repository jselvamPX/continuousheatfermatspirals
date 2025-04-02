cmake_minimum_required(VERSION 3.16)
project(Fermat3D VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 17)

# Set installation paths
include(GNUInstallDirs)
set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE PATH "Installation directory for libraries")
set(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR} CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE PATH "Installation directory for header files")

# Try to find system-installed dependencies first
find_package(Eigen3 3.4 QUIET)
if(NOT Eigen3_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
    )
    FetchContent_MakeAvailable(eigen)
endif()

# OpenMesh
find_package(OpenMesh QUIET)
if(NOT OpenMesh_FOUND)
    FetchContent_Declare(
        openmesh
        GIT_REPOSITORY https://www.graphics.rwth-aachen.de:9000/OpenMesh/OpenMesh.git
        GIT_TAG OpenMesh-11.0
        GIT_SHALLOW TRUE
    )
    set(BUILD_APPS OFF CACHE BOOL "Build OpenMesh apps")
    set(OPENMESH_DOCS OFF CACHE BOOL "Build OpenMesh documentation")
    set(BUILD_TESTING OFF CACHE BOOL "Build OpenMesh tests")
    FetchContent_MakeAvailable(openmesh)
endif()

# libigl
find_package(libigl QUIET)
if(NOT libigl_FOUND)
    FetchContent_Declare(
        libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG v2.5.0
    )
    FetchContent_MakeAvailable(libigl)
endif()

# Python bindings
if(BUILD_PYTHON_BINDINGS)
    find_package(Python3 COMPONENTS Interpreter Development Development.Module Development.Embed NumPy REQUIRED)
    if(NOT Python3_FOUND)
        message(FATAL_ERROR "Python3 not found")
    endif()
    if(NOT Python3_NumPy_FOUND)
        message(FATAL_ERROR "NumPy not found")
    endif()

    # pybind11
    find_package(pybind11 CONFIG QUIET)
    if(NOT pybind11_FOUND)
        FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11.git
            GIT_TAG v2.11.1
        )
        FetchContent_MakeAvailable(pybind11)
    endif()

    # Add Python module
    pybind11_add_module(_heat_fermat_3d_core
        heat_fermat_3d/src/heat_fermat_3d_bindings.cpp
        heat_fermat_3d/src/Spirals.cpp
        heat_fermat_3d/src/SpiralDt.cpp
    )

    target_include_directories(_heat_fermat_3d_core PRIVATE 
        ${CMAKE_CURRENT_SOURCE_DIR}/heat_fermat_3d/src
        $<$<BOOL:${OpenMesh_FOUND}>:${OpenMesh_INCLUDE_DIRS}>
        $<$<NOT:$<BOOL:${OpenMesh_FOUND}>>:${openmesh_SOURCE_DIR}/src>
        $<$<NOT:$<BOOL:${OpenMesh_FOUND}>>:${openmesh_BINARY_DIR}/src>
        ${Python3_NumPy_INCLUDE_DIRS}
    )

    target_link_libraries(_heat_fermat_3d_core PRIVATE
        igl::core
        $<$<BOOL:${OpenMesh_FOUND}>:${OpenMesh_LIBRARIES}>
        $<$<NOT:$<BOOL:${OpenMesh_FOUND}>>:OpenMeshCore>
        $<$<NOT:$<BOOL:${OpenMesh_FOUND}>>:OpenMeshTools>
        Eigen3::Eigen
    )

    # Set the output name of the Python module to ensure it's installed with the correct name
    set_target_properties(_heat_fermat_3d_core PROPERTIES
        PREFIX ""
        OUTPUT_NAME "_heat_fermat_3d_core"
        SUFFIX ".so"
    )
    
    # Install Python module and its dependencies
    install(TARGETS _heat_fermat_3d_core DESTINATION .)
    if(NOT OpenMesh_FOUND)
        # Install OpenMesh libraries
        install(TARGETS OpenMeshCore OpenMeshTools
            RUNTIME DESTINATION .
            LIBRARY DESTINATION .
            ARCHIVE DESTINATION .
        )
        # Set RPATH for the extension module
        set_target_properties(_heat_fermat_3d_core PROPERTIES
            BUILD_WITH_INSTALL_RPATH TRUE
            INSTALL_RPATH "$ORIGIN"
        )
    endif()
endif()
