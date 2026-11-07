add_subdirectory(util)
add_subdirectory(solver)
add_subdirectory(inout)

add_executable(ground_seek)
target_sources(ground_seek PRIVATE
  ground_seek.f90)

target_link_libraries(ground_seek PRIVATE
  util
  solver
  inout
  OpenMP::OpenMP_Fortran
)

