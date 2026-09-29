cmake_minimum_required(VERSION 3.23)

project(quantumDecaf
  DESCRIPTION "FOO BAR"
  LANGUAGES Fortran)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
add_subdirectory(src)

# Disable Position Independent Code globally (creates statically linked code)
set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

# Enable Link Time Optimization (LTO / IPO) for Release builds only.
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
"$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<CONFIG:Debug>>:-O0;-g;-Wall;-Wextra;-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow>"
"$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,GNU>,$<CONFIG:Release>>:-O3>"

# --- Intel (ifort / ifx) Flags ---
"$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>,$<CONFIG:Debug>>:-O0;-g;-warn;all;-check;all;-traceback;-fpe0>"
"$<$<AND:$<COMPILE_LANG_AND_ID:Fortran,Intel,IntelLLVM>,$<CONFIG:Release>>:-O3>"
)

option(BUILD_TESTING "build and run tests" ON)

include(cmake/util.cmake)

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
