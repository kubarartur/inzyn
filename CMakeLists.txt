cmake_minimum_required(VERSION 3.23)

project(quantumDecaf
  DESCRIPTION "FOO BAR"
  LANGUAGES Fortran)

add_subdirectory(src)

option(NO_CONTACT "disable contact potential" OFF)
option(NO_DIPOLE  "disable dipole-dipole potential" OFF)
option(NO_OSCILLATOR "disable external trap potential" OFF)

include(cmake/util.cmake)
