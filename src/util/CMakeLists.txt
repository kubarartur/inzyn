add_library(util)
target_sources(util PRIVATE
  phys_const.F90
  typedef.f90
  grid.f90
  initial_conditions.f90
)

set_target_properties(util PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)

target_include_directories(util INTERFACE
    ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)

