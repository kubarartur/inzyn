add_library(phys_const) target_sources(phys_const PRIVATE
  phys_const.F90
)

set_target_properties(util PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)
# so that anything importing this lib inherits
# the headers
target_include_directories(util INTERFACE
    ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)
