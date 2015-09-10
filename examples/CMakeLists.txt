include_directories("../")

add_executable(hawp_gradient hawp_gradient.cpp)
target_link_libraries(hawp_gradient yaml-cpp waveblocks-static)

add_executable(hawp_homogen hawp_homogen.cpp)

add_executable(hawp_inhomogen hawp_inhomogen.cpp)

add_executable(hawp_eval hawp_eval.cpp)