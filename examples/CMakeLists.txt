include_directories("../")

add_executable(hawp_gradient hawp_gradient.cpp)
target_link_libraries(hawp_gradient yaml-cpp waveblocks-static)

add_executable(hawp_homogen hawp_homogen.cpp)
target_link_libraries(hawp_homogen yaml-cpp waveblocks-static)

add_executable(hawp_inhomogen hawp_homogen.cpp)
target_link_libraries(hawp_homogen yaml-cpp waveblocks-static)

add_executable(hawp_performance hawp_performance.cpp)
target_link_libraries(hawp_homogen yaml-cpp waveblocks-static)