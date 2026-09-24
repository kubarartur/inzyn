add_library(solver)
target_sources(solver PRIVATE
  solver.f90
)

target_link_libraries(solver PUBLIC
  util
)

set_target_properties(solver PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)

target_include_directories(solver INTERFACE
    ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)
