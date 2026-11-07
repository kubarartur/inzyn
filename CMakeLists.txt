cmake_minimum_required(VERSION 3.23)

project(inzyn
  DESCRIPTION "FOO BAR"
  LANGUAGES Fortran)

#for autogeneration
find_package(Python3 REQUIRED COMPONENTS Interpreter)
#find fypp
find_program(fypp fypp REQUIRED)

#to satisfy fortran stdlib's requirements
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(fortran_stdlib REQUIRED)

#for -parallel ifortran flag
find_package(OpenMP REQUIRED)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE OFF)

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT error)
if(ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
else()
    message(WARNING "IPO/LTO is not supported by this compiler: ${error}")
endif()

add_compile_options(
#gfortran
"$<$<AND:$<Fortran_COMPILER_ID:GNU>,$<CONFIG:Debug>>:-O0;-g;-Wall;-Wextra;-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow>"
"$<$<AND:$<Fortran_COMPILER_ID:GNU>,$<CONFIG:Release>>:-O3>"

#ifort/ifx
"$<$<AND:$<Fortran_COMPILER_ID:Intel,IntelLLVM>,$<CONFIG:Release>>:-parallel>"
)

add_subdirectory(src)

option(BUILD_TESTING "build and run tests" ON)

include(cmake/util.cmake)

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()
