
add_executable (vpic-reader vpic-reader.cxx)

foreach (lcv ${MPI_CXX_COMPILE_FLAGS_LIST})
    if (NOT ${lcv} STREQUAL "")
        target_compile_options (vpic-reader PUBLIC ${lcv})
    endif ()
endforeach ()

# XXX: have to do this one dir at a time otherwise, otherwise I get
# error: target 'deltafs' INTERFACE_INCLUDE_DIRECTORIES contains path
#               prefixed in the source directory
foreach (lcv ${MPI_CXX_INCLUDE_PATH})
    target_include_directories (vpic-reader PUBLIC ${lcv})
endforeach ()

foreach (lcv ${MPI_CXX_LIBRARIES})
    target_link_libraries (vpic-reader ${lcv})
endforeach ()

set_property (TARGET vpic-reader APPEND PROPERTY LINK_FLAGS
              ${MPI_CXX_LINK_FLAGS})
#
# "make install" rules
#
install (TARGETS vpic-reader
         RUNTIME DESTINATION bin)
