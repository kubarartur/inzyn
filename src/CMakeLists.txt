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

add_executable(real_time_evol)
target_sources(real_time_evol PRIVATE
  real_time_evol.f90)

target_link_libraries(real_time_evol PRIVATE
  util
  solver
  inout
  OpenMP::OpenMP_Fortran
)
