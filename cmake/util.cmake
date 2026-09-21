function(add_logging_test testname fname)
    add_executable(${testname})
    target_sources(${testname} 
        PRIVATE
            ${fname}
    )

    add_test(
        NAME ${testname}
        COMMAND $<TARGET_FILE:${testname}>
    )
endfunction()
