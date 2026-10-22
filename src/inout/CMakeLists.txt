add_library(inout)
target_sources(inout PRIVATE
  initializer.f90
  arr_iout.f90
)
target_link_libraries(inout PUBLIC
  util
  solver
)

set_target_properties(inout PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)

target_include_directories(inout INTERFACE
  ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)
