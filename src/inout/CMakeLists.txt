add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/arr_iout.f90
  COMMAND ${fypp} ${CMAKE_CURRENT_SOURCE_DIR}/arr_iout.fypp ${CMAKE_CURRENT_BINARY_DIR}/arr_iout.f90
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/arr_iout.fypp
  COMMENT "Calling preprocessor on arr_inout.fypp"
)

add_library(inout)
target_sources(inout PRIVATE
  initializer.f90
  ${CMAKE_CURRENT_BINARY_DIR}/arr_iout.f90
  functor/double_gaussian.f90
  functor/double_well.f90
  functor/single_well.f90
  functor/double_well_shift.f90
  functor/single_well_shift.f90
  functor/from_file.f90
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
