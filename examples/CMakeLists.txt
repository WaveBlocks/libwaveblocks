include_directories("../")

#add_executable(hawp_gradient hawp_gradient.cpp)
#target_link_libraries(hawp_gradient yaml-cpp waveblocks-static)

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
