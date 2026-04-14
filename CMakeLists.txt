cmake_minimum_required(VERSION 3.16)
project(power-knitting LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Prevent in-source builds
if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same.")
endif()

# Output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ── Compiler flags ────────────────────────────────────────────────────────────
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(DISABLED_WARNINGS
        "-Wno-unused-parameter -Wno-unused-variable -Wno-unused-function \
         -Wno-deprecated-declarations -Wno-missing-braces")

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "-Wall -g3 -ferror-limit=5 -fcolor-diagnostics ${DISABLED_WARNINGS}")
        set(CMAKE_CXX_FLAGS_DEBUG "-fsanitize=address -fno-limit-debug-info -g3")
    else()
        set(CMAKE_CXX_FLAGS "-Wall -g3 -fmax-errors=5 ${DISABLED_WARNINGS} -Wno-maybe-uninitialized")
        set(CMAKE_CXX_FLAGS_DEBUG "-g3")
    endif()

    set(CMAKE_CXX_FLAGS_RELEASE        "-march=native -O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
endif()

# ── Dependencies via FetchContent ─────────────────────────────────────────────
include(FetchContent)

# Eigen3: system preferred, FetchContent fallback
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found via find_package — fetching 3.4.0")
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        5.0
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(eigen)
endif()

# CLI11
FetchContent_Declare(
    CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG        v2.6.2
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(CLI11)

# geometry-central (custom fork — submodule at deps/geometry-central)
add_subdirectory(deps/geometry-central)

# nlohmann/json — declared before polyscope so polyscope reuses this target
FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.3
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(nlohmann_json)

# Polyscope
FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG        v2.6.1
)
FetchContent_MakeAvailable(polyscope)

# libigl (header-only core)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG        v2.5.0
    GIT_SHALLOW    TRUE
)
set(LIBIGL_USE_STATIC_LIBRARY OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(libigl)

# ── OpenMP (Apple Silicon / Homebrew libomp) ──────────────────────────────────
if(APPLE)
    set(OpenMP_C_FLAGS        "-Xpreprocessor -fopenmp")
    set(OpenMP_CXX_FLAGS      "-Xpreprocessor -fopenmp")
    set(OpenMP_C_LIB_NAMES    "omp")
    set(OpenMP_CXX_LIB_NAMES  "omp")
    set(OpenMP_omp_LIBRARY    omp)
    include_directories("/opt/homebrew/opt/libomp/include")
    link_directories("/opt/homebrew/opt/libomp/lib")
endif()

find_package(OpenMP)
if(OpenMP_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(STATUS "OpenMP NOT found — building without parallelization")
endif()

# ── Executable ────────────────────────────────────────────────────────────────
add_executable(power_knitting src/main.cpp src/knit_model.cpp src/time_function.cpp)

target_link_libraries(power_knitting PRIVATE
    polyscope
    geometry-central
    Eigen3::Eigen
    CLI11::CLI11
    nlohmann_json::nlohmann_json
    igl::core
    ${OpenMP_CXX_LIBRARIES}
)
