include_directories("../")

# add_executable(hawp_gradient hawp_gradient.cpp)
# target_link_libraries(hawp_gradient yaml-cpp waveblocks-static)

add_executable(hawp_homogen hawp_homogen.cpp)

add_executable(hawp_inhomogen hawp_inhomogen.cpp)

add_executable(hawp_eval hawp_eval.cpp)

#add_executable(hawp_t hawp_t.cpp)

add_executable(anharmonic_1D anharmonic_1D.cpp)
target_link_libraries(anharmonic_1D ${LINK_LIBS})

add_executable(harmonic_2D harmonic_2D.cpp)
target_link_libraries(harmonic_2D ${LINK_LIBS})

add_executable(tunneling_1D tunneling_1D.cpp)
target_link_libraries(tunneling_1D ${LINK_LIBS})

add_executable(tunneling2_1D tunneling2_1D.cpp)
target_link_libraries(tunneling2_1D ${LINK_LIBS})

add_executable(nonsense_2N_1D nonsense_2N_1D.cpp)
add_executable(nonsense_2N_1D_inhom nonsense_2N_1D_inhom.cpp)
add_executable(nonsense_2N_2D nonsense_2N_2D.cpp)
add_executable(nonsense_2N_2D_inhom nonsense_2N_2D_inhom.cpp)


############################################
# Testing new propagators
############################################

add_executable(anharmonic_1D_newpropagator anharmonic_1D_newpropagator.cpp)
target_link_libraries(anharmonic_1D_newpropagator ${LINK_LIBS})

add_executable(harmonic_2D_newpropagator harmonic_2D_newpropagator.cpp)
target_link_libraries(harmonic_2D_newpropagator ${LINK_LIBS})

add_executable(tunneling_1D_newpropagator tunneling_1D_newpropagator.cpp)
target_link_libraries(tunneling_1D_newpropagator ${LINK_LIBS})

add_executable(tunneling2_1D_newpropagator tunneling2_1D_newpropagator.cpp)
target_link_libraries(tunneling2_1D_newpropagator ${LINK_LIBS})

add_executable(nonsense_2N_1D_newpropagator nonsense_2N_1D_newpropagator.cpp)
target_link_libraries(nonsense_2N_1D_newpropagator ${LINK_LIBS})

add_executable(nonsense_2N_1D_inhom_newpropagator nonsense_2N_1D_inhom_newpropagator.cpp)
target_link_libraries(nonsense_2N_1D_inhom_newpropagator ${LINK_LIBS})

add_executable(nonsense_2N_2D_newpropagator nonsense_2N_2D_newpropagator.cpp)
target_link_libraries(nonsense_2N_2D_newpropagator ${LINK_LIBS})

add_executable(nonsense_2N_2D_inhom_newpropagator nonsense_2N_2D_inhom_newpropagator.cpp)
target_link_libraries(nonsense_2N_2D_inhom_newpropagator ${LINK_LIBS})

add_executable(harmonic_XD_newpropagator harmonic_XD_newpropagator.cpp)
target_link_libraries(harmonic_XD_newpropagator ${LINK_LIBS})

add_executable(henon-heiles_2D_newpropagator henon-heiles_2D_newpropagator.cpp)
target_link_libraries(henon-heiles_2D_newpropagator ${LINK_LIBS})
