file(GLOB reg_targets_sources *.f90)

foreach(target_source ${reg_targets_sources})
    get_filename_component( bfname "${target_source}" NAME_WLE)
    message("test name is " ${bfname})
    add_logging_test(${bfname} ${target_source})
    target_link_libraries( ${bfname}
        PRIVATE
            util
    )
endforeach()

file(GLOB fpp_targets_sources *.F90)

foreach(target_source ${fpp_targets_sources})
    get_filename_component( bfname "${target_source}" NAME_WLE)
    message("test name is " ${bfname})
    add_logging_test(${bfname} ${target_source})
    target_link_libraries( ${bfname}
        PRIVATE
            util
    )
    target_compile_definitions(${bfname}
        PRIVATE
            dumpname="${CMAKE_BINARY_DIR}/position_dump/${bfname}.dat")
endforeach()
