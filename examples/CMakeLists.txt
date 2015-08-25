include_directories("../")

add_executable(hawp_homogen_hard hawp_homogen_hard.cpp)
add_executable(hawp_homogen_easy hawp_homogen_easy.cpp)

add_executable(hawp_gradient_easy hawp_gradient_easy.cpp)
target_link_libraries(hawp_gradient_easy yaml-cpp waveblocks-static)