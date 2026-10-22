add_library(testmodules)
target_sources(testmodules PRIVATE
  mod_test_lap.f90
)

target_link_libraries(testmodules PUBLIC
  util
)
set_target_properties(testmodules PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)

target_include_directories(testmodules INTERFACE
    ${CMAKE_CURRENT_BINARY_DIR}/mod_files
)

