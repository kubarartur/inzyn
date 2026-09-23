cmake_minimum_required(VERSION 3.23)

project(quantumDecaf
  DESCRIPTION "FOO BAR"
  LANGUAGES Fortran)

add_subdirectory(src)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(
        CMAKE_Fortran_FLAGS_INIT
        "-fimplicit-none"
        "-ffree-line-length-132"
    )
    set(
        CMAKE_Fortran_FLAGS_RELEASE_INIT
    )


    block()
        #contents of wall, without Wunused and Winteger-division
        #TODO: consider chaining to -Wall -Wno-unused Wno-integer-division
        set(wall_ish "-fsanitize=address;-Waliasing;-Wampersand;-Wconversion;-Wsurprising;-Wc-binding-type;-Wintrinsics-std;-Wtabs;-Wintrinsic-shadow;-Wline-truncation;-Wtarget-lifetime;-Wreal-q-constant;-Wundefined-do-loop;-Wmaybe-uninitialized")
        set(loc_deb_flags "-O0;-g;-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow;-finit-real=snan;-Werror;-pedantic;-fimplicit-none;-std=f2018")
        list(APPEND loc_deb_flags ${wall_ish})
        set(locl_deb_flags "-fsanitize=address;-O0;-g;-fcheck=all;-fbacktrace;-ffpe-trap=invalid,zero,overflow;-finit-real=snan")
        set(loclr_el_flags "-O3;-march=native;-flto;-funroll-loops")

        set(loc_rel_flags "-O3;-march=native;-flto;-funroll-loops;-Wall;-Wextra;-pedantic;-fimplicit-none;-std=f2018")
        foreach(flag ${loc_deb_flags})
            add_compile_options(
                $<$<CONFIG:Debug>:${flag}>
            )
        endforeach()

        foreach(flag ${locl_deb_flags})
            add_link_options(
                $<$<CONFIG:Debug>:${flag}>
            )
        endforeach()

        foreach(flag ${loc_rel_flags})
            add_compile_options(
                $<$<CONFIG:Release>:${flag}>
            )
        endforeach()

        foreach(flag ${locl_rel_flags})
            add_link_options(
                $<$<CONFIG:Release>:${flag}>
            )
        endforeach()
    endblock()

else()
    message(WARNING "can't setup proper compilation flags for used compiler.
                     defaulting to basic flags")
endif()

option(NO_CONTACT "disable contact potential" OFF)
option(NO_DIPOLE  "disable dipole-dipole potential" OFF)
option(NO_OSCILLATOR "disable external trap potential" OFF)

include(cmake/util.cmake)
