add_library(solver)
target_sources(solver PRIVATE
  solver.f90
  hamiltonian.f90
  autogen_laplaceian.f90
  autogen_vdd_boundary.f90
)
find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_custom_command(
  OUTPUT autogen_laplaceian.f90
  COMMAND Python3::Interpreter "${CMAKE_CURRENT_SOURCE_DIR}/laplgen.py" autogen_laplaceian.f90
  DEPENDS laplgen.py
  COMMENT "Generating submodule with Laplaceian Calculation"
  VERBATIM
)

add_custom_command(
  OUTPUT autogen_vdd_boundary.f90
  COMMAND Python3::Interpreter "${CMAKE_CURRENT_SOURCE_DIR}/vdd_intgen.py" autogen_vdd_boundary.f90
  DEPENDS vdd_intgen.py
  COMMENT "Generating submodule with Laplaceian Calculation"
  VERBATIM
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
