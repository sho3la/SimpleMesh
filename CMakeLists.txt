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

# The header-only exact-arithmetic / arrangement kernel. Header-only, so these
# are attached to the target purely so they appear in IDE project trees.
set(SIMPLEMESH_EXACT_HEADERS
    include/simplemesh/exact/Expansion.h
    include/simplemesh/exact/Interval.h
    include/simplemesh/exact/Predicates.h
    include/simplemesh/exact/HomogeneousGeometry.h
    include/simplemesh/exact/CDT2d.h
    include/simplemesh/exact/TriangleIntersection.h
    include/simplemesh/exact/MeshInTriangle.h
    include/simplemesh/exact/SurfaceIntersection.h
    include/simplemesh/exact/MeshBoolean.h)

add_library(simplemesh_core STATIC
    ${SIMPLEMESH_SOURCES}
    ${SIMPLEMESH_HEADERS}
    ${SIMPLEMESH_EXACT_HEADERS})

target_include_directories(simplemesh_core
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Mirror the on-disk layout in IDE project trees ("Source Files" / "Header
# Files" folders instead of one flat list).
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/src
    PREFIX "Source Files" FILES ${SIMPLEMESH_SOURCES})
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/include
    PREFIX "Header Files" FILES ${SIMPLEMESH_HEADERS} ${SIMPLEMESH_EXACT_HEADERS})

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

    # Mesh fixtures read by the tests live under the top-level data/ folder.
    # Round-trip I/O tests write their scratch meshes into a build-tree scratch
    # dir (SM_SCRATCH_DIR), never the source tree / repo root.
    set(SM_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
    set(SM_SCRATCH_DIR "${CMAKE_BINARY_DIR}/test_scratch")
    file(MAKE_DIRECTORY "${SM_SCRATCH_DIR}")

    add_executable(simplemesh_tests tests/test_mesh.cpp)
    target_link_libraries(simplemesh_tests PRIVATE simplemesh_core)
    target_compile_definitions(simplemesh_tests PRIVATE
        SM_SCRATCH_DIR="${SM_SCRATCH_DIR}")
    add_test(NAME simplemesh_tests COMMAND simplemesh_tests
             WORKING_DIRECTORY "${SM_SCRATCH_DIR}")

    # --- exact self-intersection / arrangement: layered tests ---------------
    # Baseline oracle: validates the test corpus and records the "before"
    # self-intersection counts the end-to-end stage must drive to zero.
    add_executable(test_corpus_baseline tests/exact/test_corpus_baseline.cpp)
    target_link_libraries(test_corpus_baseline PRIVATE simplemesh_core)
    target_compile_definitions(test_corpus_baseline PRIVATE
        SM_ISECT_DATA_DIR="${SM_DATA_DIR}/isect")
    add_test(NAME test_corpus_baseline COMMAND test_corpus_baseline)

    # Exact arithmetic must NOT be contracted into FMAs, or the error-free
    # two_sum/two_product transforms break. Apply to every exact-kernel test.
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(SM_EXACT_FLAGS -ffp-contract=off)
    elseif(MSVC)
        set(SM_EXACT_FLAGS /fp:strict)
    else()
        set(SM_EXACT_FLAGS "")
    endif()

    # Multi-precision Expansion. Header-only; validated against an
    # exact 128-bit integer oracle.
    add_executable(test_expansion tests/exact/test_expansion.cpp tests/exact/Int128.h)
    target_link_libraries(test_expansion PRIVATE simplemesh_core)
    target_compile_options(test_expansion PRIVATE ${SM_EXACT_FLAGS})
    add_test(NAME test_expansion COMMAND test_expansion)

    # Round-to-nearest interval filter. Validated for the filter safety
    # invariant (a determined sign is never wrong) against the exact 128-bit integer oracle.
    add_executable(test_interval tests/exact/test_interval.cpp tests/exact/Int128.h)
    target_link_libraries(test_interval PRIVATE simplemesh_core)
    target_compile_options(test_interval PRIVATE ${SM_EXACT_FLAGS})
    add_test(NAME test_interval COMMAND test_interval)

    # Exact predicates orient2d / orient3d / in_circle. Golden-tested
    # against an exact 128-bit determinant oracle (zero tolerance for mismatch).
    add_executable(test_predicates tests/exact/test_predicates.cpp tests/exact/Int128.h)
    target_link_libraries(test_predicates PRIVATE simplemesh_core)
    target_compile_options(test_predicates PRIVATE ${SM_EXACT_FLAGS})
    add_test(NAME test_predicates COMMAND test_predicates)

    # Exact homogeneous geometry (vec2HE/vec3HE) + projective predicates.
    add_executable(test_homogeneous tests/exact/test_homogeneous.cpp tests/exact/Int128.h)
    target_link_libraries(test_homogeneous PRIVATE simplemesh_core)
    target_compile_options(test_homogeneous PRIVATE ${SM_EXACT_FLAGS})
    add_test(NAME test_homogeneous COMMAND test_homogeneous)

    # Constrained Delaunay triangulation - combinatorial core + point
    # insertion. Validated by the exact empty-circumcircle Delaunay property.
    add_executable(test_cdt tests/exact/test_cdt.cpp)
    target_link_libraries(test_cdt PRIVATE simplemesh_core)
    target_compile_options(test_cdt PRIVATE ${SM_EXACT_FLAGS})
    add_test(NAME test_cdt COMMAND test_cdt)

    # CDT constraint insertion (recovery + crossing intersections).
    add_executable(test_cdt_constraints tests/exact/test_cdt_constraints.cpp)
    target_link_libraries(test_cdt_constraints PRIVATE simplemesh_core)
    target_compile_options(test_cdt_constraints PRIVATE ${SM_EXACT_FLAGS})
    add_test(NAME test_cdt_constraints COMMAND test_cdt_constraints)

    # Symbolic triangle/triangle intersection (TriangleRegion).
    add_executable(test_triangle_isect tests/exact/test_triangle_isect.cpp)
    target_link_libraries(test_triangle_isect PRIVATE simplemesh_core)
    target_compile_options(test_triangle_isect PRIVATE ${SM_EXACT_FLAGS})
    add_test(NAME test_triangle_isect COMMAND test_triangle_isect)

    # ExactCDT2d (exact-coordinate CDT) + MeshInTriangle per-facet remesh.
    add_executable(test_mesh_in_triangle tests/exact/test_mesh_in_triangle.cpp)
    target_link_libraries(test_mesh_in_triangle PRIVATE simplemesh_core)
    target_compile_options(test_mesh_in_triangle PRIVATE ${SM_EXACT_FLAGS})
    add_test(NAME test_mesh_in_triangle COMMAND test_mesh_in_triangle)

    # Arrangement assembly - broad phase + symbolic intersection +
    # per-facet exact remesh + exact global vertex merge. Verifies on the corpus
    # that ZERO non-adjacent triangle pairs still interpenetrate (exact oracle).
    add_executable(test_surface_intersection tests/exact/test_surface_intersection.cpp)
    target_link_libraries(test_surface_intersection PRIVATE simplemesh_core)
    target_compile_options(test_surface_intersection PRIVATE ${SM_EXACT_FLAGS})
    target_compile_definitions(test_surface_intersection PRIVATE
        SM_ISECT_DATA_DIR="${SM_DATA_DIR}/isect")
    add_test(NAME test_surface_intersection COMMAND test_surface_intersection)

    # Shell classification & boolean ops (union/intersection/difference) via the
    # Weiler model + radial sort. Validated on cube_x_cube and nested cubes by
    # watertight + manifold + zero self-intersection + exact volume + bounding box.
    add_executable(test_mesh_boolean tests/exact/test_mesh_boolean.cpp)
    target_link_libraries(test_mesh_boolean PRIVATE simplemesh_core)
    target_compile_options(test_mesh_boolean PRIVATE ${SM_EXACT_FLAGS})
    target_compile_definitions(test_mesh_boolean PRIVATE
        SM_ISECT_DATA_DIR="${SM_DATA_DIR}/isect")
    add_test(NAME test_mesh_boolean COMMAND test_mesh_boolean)
endif()
