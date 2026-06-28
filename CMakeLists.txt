# ============================================================================
#  SimpleMesh - top-level CMake build script
# ----------------------------------------------------------------------------
#  Targets:
#    * simplemesh_core   : static C++ library (the data structure + algorithms)
#    * simplemesh        : the Python extension module (pybind11)
#    * cube_example      : a standalone C++ demo            (SIMPLEMESH_BUILD_APPS)
#    * simplemesh_tests  : the C++ self-test                (SIMPLEMESH_BUILD_APPS)
#
#  Build standalone:
#      cmake -S . -B build
#      cmake --build build --config Release
#  Or install the Python package (uses pyproject.toml + scikit-build-core):
#      pip install .
# ============================================================================
cmake_minimum_required(VERSION 3.18)

project(SimpleMesh
    VERSION 0.1.0
    DESCRIPTION "Simplified halfedge mesh library"
    LANGUAGES CXX)

# ----- global settings ------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(SIMPLEMESH_BUILD_PYTHON "Build the Python bindings"        ON)
option(SIMPLEMESH_BUILD_APPS   "Build the demo + test executables" ON)

# ----------------------------------------------------------------------------
#  1. The core C++ library
# ----------------------------------------------------------------------------
set(SIMPLEMESH_SOURCES
    src/Mesh.cpp
    src/MeshIO.cpp
    src/Algorithms.cpp
    src/Decimater.cpp
    src/MeshChecker.cpp
    src/MeshRepair.cpp
    src/PreMesh.cpp)

# Listed in the target so they show up in IDE project trees (Visual Studio
# Solution Explorer, Xcode, etc.). They are headers, so nothing is compiled
# from them directly.
set(SIMPLEMESH_HEADERS
    include/simplemesh/Vec3.h
    include/simplemesh/Handles.h
    include/simplemesh/Items.h
    include/simplemesh/Properties.h
    include/simplemesh/Mesh.h
    include/simplemesh/Circulators.h
    include/simplemesh/SmartHandles.h
    include/simplemesh/Algorithms.h
    include/simplemesh/Decimater.h
    include/simplemesh/MeshChecker.h
    include/simplemesh/MeshRepair.h
    include/simplemesh/PreMesh.h)

add_library(simplemesh_core STATIC
    ${SIMPLEMESH_SOURCES}
    ${SIMPLEMESH_HEADERS})

target_include_directories(simplemesh_core
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Mirror the on-disk layout in IDE project trees ("Source Files" / "Header
# Files" folders instead of one flat list).
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/src
    PREFIX "Source Files" FILES ${SIMPLEMESH_SOURCES})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/include
    PREFIX "Header Files" FILES ${SIMPLEMESH_HEADERS})

# Position-independent code so the static lib links into the shared Python
# module on Linux/Mac. Harmless on Windows.
set_target_properties(simplemesh_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ----------------------------------------------------------------------------
#  2. Python bindings (pybind11). Prefer an installed pybind11 (e.g. provided by
#     scikit-build-core during `pip install`); otherwise fetch it.
# ----------------------------------------------------------------------------
if(SIMPLEMESH_BUILD_PYTHON)
    find_package(pybind11 CONFIG QUIET)
    if(NOT pybind11_FOUND)
        include(FetchContent)
        FetchContent_Declare(
            pybind11
            GIT_REPOSITORY https://github.com/pybind/pybind11.git
            GIT_TAG        v2.13.6)
        FetchContent_MakeAvailable(pybind11)
    endif()

    pybind11_add_module(simplemesh python/bindings.cpp)
    target_link_libraries(simplemesh PRIVATE simplemesh_core)

    # For `pip install` (scikit-build-core): place the module at the wheel root.
    install(TARGETS simplemesh LIBRARY DESTINATION . RUNTIME DESTINATION .)
endif()

# ----------------------------------------------------------------------------
#  3. C++ examples + self-test (skipped for wheel builds)
# ----------------------------------------------------------------------------
if(SIMPLEMESH_BUILD_APPS)
    # One executable per examples/NN_*.cpp (mirrors the Python examples).
    set(SIMPLEMESH_EXAMPLES
        01_build_mesh  02_io          03_navigation  04_smart_handles
        05_editing     06_properties  07_buffers     08_geometry
        09_subdivision 10_decimation  11_smoothing   12_check_and_repair)
    foreach(ex ${SIMPLEMESH_EXAMPLES})
        add_executable(${ex} examples/${ex}.cpp)
        target_link_libraries(${ex} PRIVATE simplemesh_core)
    endforeach()

    enable_testing()
    add_executable(simplemesh_tests tests/test_mesh.cpp)
    target_link_libraries(simplemesh_tests PRIVATE simplemesh_core)
    add_test(NAME simplemesh_tests COMMAND simplemesh_tests)
endif()
