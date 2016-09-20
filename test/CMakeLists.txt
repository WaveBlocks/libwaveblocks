include_directories("../")

add_executable(test_eigen test_eigen.cpp)
target_link_libraries(test_eigen ${LINK_LIBS})

#add_executable(test_eigenhdf test_eigenhdf.cpp)
#target_link_libraries(test_eigenhdf ${LINK_LIBS})

add_executable(test_yaml test_yaml.cpp)
target_link_libraries(test_yaml waveblocks-static ${LINK_LIBS})

add_executable(test_yaml_eigen test_yaml_eigen.cpp)
target_link_libraries(test_yaml_eigen waveblocks-static ${LINK_LIBS})

add_executable(test_yaml_pi test_yaml_pi.cpp)
target_link_libraries(test_yaml_pi waveblocks-static ${LINK_LIBS})

#add_executable(test_yaml_hdf5 test_yaml_hdf5.cpp)
#target_link_libraries(test_yaml_hdf5 waveblocks-static ${LINK_LIBS})

add_executable(test_cont_sqrt test_cont_sqrt.cpp)

add_executable(test_param_mixing test_param_mixing.cpp)

add_executable(test_innerproduct test_innerproduct.cpp)
add_executable(test_innerproduct_runtimes test_innerproduct_runtimes.cpp)
add_executable(test_innerproduct_problems test_innerproduct_problems.cpp)

add_executable(test_matexp test_matexp.cpp)

if(GTEST_FOUND)
    add_executable(datatest datatest.cpp)
    target_link_libraries(datatest ${LINK_LIBS})
endif()
