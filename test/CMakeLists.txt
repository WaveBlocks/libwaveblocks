include_directories("../")

add_executable(test_cont_sqrt test_cont_sqrt.cpp)

add_executable(test_eigen test_eigen.cpp)
target_link_libraries(test_eigen ${LINK_LIBS})

add_executable(test_eigenhdf test_eigenhdf.cpp)
target_link_libraries(test_eigenhdf ${LINK_LIBS})

add_executable(test_yaml test_yaml.cpp)
target_link_libraries(test_yaml waveblocks-static ${LINK_LIBS})

add_executable(test_yaml_eigen test_yaml_eigen.cpp)
target_link_libraries(test_yaml_eigen waveblocks-static ${LINK_LIBS})

add_executable(test_innerproduct test_innerproduct.cpp)
