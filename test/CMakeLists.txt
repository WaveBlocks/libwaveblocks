include_directories("../")

add_executable(test_cont_sqrt test_cont_sqrt.cpp)

add_executable(test_innerproduct
    test_innerproduct.cpp
)

add_executable(test_eigenhdf test_eigenhdf.cpp)
target_link_libraries(test_eigenhdf ${LINK_LIBS})
