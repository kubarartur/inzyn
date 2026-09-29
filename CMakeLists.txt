cmake_minimum_required(VERSION 3.23)

project(quantumDecaf
  DESCRIPTION "FOO BAR"
  LANGUAGES Fortran)

add_subdirectory(src)

# Disable Position Independent Code globally (creates statically linked code)
set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

# Enable Link Time Optimization (LTO / IPO) for Release builds only.
# CMake natively abstracts the compiler flag (e.g., -flto for GNU, -ipo for Intel).
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)

# Optional but recommended: Check if the current compiler actually supports IPO
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT error)
if(ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
else()
    message(WARNING "IPO/LTO is not supported by this compiler: ${error}")
endif()

# ---------------------------------------------------------------------------
# 2. Compiler-Specific Flags via Generator Expressions
# ---------------------------------------------------------------------------

add_compile_options(
    # --- GNU (gfortran) Flags ---
    "$<$<AND:$<Fortran_COMPILER_ID:GNU>,$<CONFIG:Debug>>:-O0;-g;-Wall;-Wextra;-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow>"
    "$<$<AND:$<Fortran_COMPILER_ID:GNU>,$<CONFIG:Release>>:-O3>"

    # --- Intel (ifort / ifx) Flags ---
    "$<$<AND:$<IN_LIST:$<Fortran_COMPILER_ID>,Intel;IntelLLVM>,$<CONFIG:Debug>>:-O0;-g;-warn;all;-check;all;-traceback;-fpe0>"
    "$<$<AND:$<IN_LIST:$<Fortran_COMPILER_ID>,Intel;IntelLLVM>,$<CONFIG:Release>>:-O3>"
)

option(BUILD_TESTING "build and run tests" ON)

include(cmake/util.cmake)

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
