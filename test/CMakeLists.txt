include_directories("../")

add_executable(test_cont_sqrt test_cont_sqrt.cpp)

add_executable(test_eigen test_eigen.cpp)
target_link_libraries(test_eigen ${LINK_LIBS})

add_executable(test_eigenhdf test_eigenhdf.cpp)
target_link_libraries(test_eigenhdf ${LINK_LIBS})

add_executable(test_param_mixing test_param_mixing.cpp)
