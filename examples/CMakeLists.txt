include_directories("../")

#add_executable(hawp_gradient hawp_gradient.cpp)
#target_link_libraries(hawp_gradient yaml-cpp waveblocks-static)

# add_executable(hawp_homogen hawp_homogen.cpp)

# add_executable(hawp_inhomogen hawp_inhomogen.cpp)

# add_executable(hawp_eval hawp_eval.cpp)

# #add_executable(hawp_t hawp_t.cpp)

# add_executable(anharmonic_1D anharmonic_1D.cpp)
# target_link_libraries(anharmonic_1D ${LINK_LIBS})

# add_executable(harmonic_2D harmonic_2D.cpp)
# target_link_libraries(harmonic_2D ${LINK_LIBS})

# add_executable(tunneling_1D tunneling_1D.cpp)
# target_link_libraries(tunneling_1D ${LINK_LIBS})

# add_executable(tunneling2_1D tunneling2_1D.cpp)
# target_link_libraries(tunneling2_1D ${LINK_LIBS})

# add_executable(nonsense_2N_1D nonsense_2N_1D.cpp)
# add_executable(nonsense_2N_1D_inhom nonsense_2N_1D_inhom.cpp)
# add_executable(nonsense_2N_2D nonsense_2N_2D.cpp)
# add_executable(nonsense_2N_2D_inhom nonsense_2N_2D_inhom.cpp)

############################################
# Testing new propagators
############################################
# TODO: remove
add_executable(anharmonic_1D_newpropagator anharmonic_1D_newpropagator.cpp)
target_link_libraries(anharmonic_1D_newpropagator ${LINK_LIBS})

# add all executables ending with _newpropagator.cpp
# file( GLOB APP_SOURCES RELATIVE app/*.cxx )
file( GLOB SRC_NEWPROPAGATOR RELATIVE "*_newpropagator.cpp" )
foreach( newexamplesrc ${SRC_NEWPROPAGATOR} )

    string( REPLACE "_newpropagator.cpp" ".cpp" oldexamplesrc ${newexamplesrc} ) # old propagator file
    string( REPLACE ".cpp" "" newexampleexe ${newexamplesrc} ) # new propagator file
    string( REPLACE ".cpp" "" oldexampleexe ${oldexamplesrc} ) # new propagator file

    add_executable( ${newexampleexe} ${newexamplesrc} )
	target_link_libraries( ${newexampleexe} ${LINK_LIBS} )

    add_executable( ${oldexampleexe} ${oldexamplesrc} )
	target_link_libraries( ${oldexampleexe} ${LINK_LIBS} )

endforeach( newexamplesrc ${SRC_NEWPROPAGATOR} )
