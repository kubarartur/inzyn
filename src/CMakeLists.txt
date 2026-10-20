add_subdirectory(util)
add_subdirectory(solver)

add_executable(main)
target_sources(main PRIVATE
  main.f90)

target_link_libraries(main PRIVATE
util
solver
OpenMP::OpenMP_Fortran
)

