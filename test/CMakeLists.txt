include_directories("../")

add_executable(slim_recusion slim_recursion.cpp util/time.cpp)
add_executable(shape_gradient shape_gradient.cpp util/time.cpp)

target_link_libraries(slim_recusion pthread)