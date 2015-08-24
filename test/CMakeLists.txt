include_directories("../")


add_executable(test_hawp_evaluator 
    test_hawp_evaluator.cpp 
    coefficients_file_parser.cpp
    complexnumber_parser.cpp
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

add_executable(yaml_test
    complexnumber_parser.cpp
    yaml_test.cpp
)

target_link_libraries(yaml_test yaml-cpp)

# add_executable(spielwiese
#     spielwiese.cpp util/time.cpp
# )