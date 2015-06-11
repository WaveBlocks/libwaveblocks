include_directories("../")


add_executable(test_hawp_evaluator 
    test_hawp_evaluator.cpp 
    util/time.cpp
)

add_executable(test_hawp_gradient 
    test_hawp_gradient.cpp 
    util/time.cpp
)

add_executable(test_shape_union
    test_shape_union.cpp
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