include_directories("../")

add_executable(test_hagedorn_wavepacket_evaluation 
    test_hagedorn_wavepacket_evaluation.cpp 
    util/time.cpp
)

add_executable(test_hagedorn_gradient_operator 
    test_hagedorn_gradient_operator.cpp 
    util/time.cpp
)