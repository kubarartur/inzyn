add_library(util)
target_sources(util PRIVATE
  typedef.f90
  phys_const.F90
  grid.f90
  system.f90
  initial_conditions.f90
  io.f90
)

set_target_properties(util PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)

target_include_directories(util INTERFACE
    ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)

