include_directories("../")


add_executable(test_hagedorn_wavepacket_evaluation 
    test_hagedorn_wavepacket_evaluation.cpp 
    util/time.cpp
)

add_executable(profile_hagedorn_wavepacket_evaluation 
    profile_hagedorn_wavepacket_evaluation.cpp 
    util/time.cpp
)

add_executable(test_hagedorn_gradient_operator 
    test_hagedorn_gradient_operator.cpp 
    util/time.cpp
)

add_executable(test_homogeneous_wavepacket
    test_homogeneous_wavepacket.cpp
    util/time.cpp
)

add_executable(test_shape_extension
    test_shape_extension.cpp
)

add_executable(test_cont_sqrt
    test_cont_sqrt.cpp
)

add_executable(spielwiese
    spielwiese.cpp util/time.cpp
)