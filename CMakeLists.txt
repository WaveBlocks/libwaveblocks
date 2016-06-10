cmake_minimum_required(VERSION 2.8)

project(waveblocks)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake_modules")

# GCC
add_definitions(
  -std=c++11
  -Wall
  -Wno-deprecated-declarations
  -Wextra
  -pedantic
  -Ofast # 03 and some optimization not valid for all standard-compliant programs (such as -ffast-math
  -fprofile-use # conditions optimizations on a profiling run
  -static # statically links library calls
  -ftree-loop-if-convert-stores # rewrites conditional writes in loops to unconiditional ones
  -ftree-loop-distribution # distributes loops for better cache usage
  )

# Clang
# add_definitions(
#  -std=c++11
#  -Wall
#  -Wno-deprecated-declarations
#  -Wextra
#  -pedantic
#  -Ofast # 03 and some optimization not valid for all standard-compliant programs (such as -ffast-math
#  -static # statically links library calls
#  -fopenmp
#  )

#add all headers to source tree to let Qtcreator know all files
file(GLOB_RECURSE WAVBLOCK_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/waveblocks
    ${CMAKE_CURRENT_SOURCE_DIR}/waveblocks/*.hpp)

add_library(waveblocks-static STATIC
  waveblocks/csv/coefficients_file_parser.cpp
  waveblocks/utilities/complexnumber_parser.cpp
  ${WAVBLOCK_HEADERS})


include(FindPkgConfig)
if(PKG_CONFIG_FOUND)
  pkg_search_module(EIGEN REQUIRED eigen3)
  if(EIGEN_FOUND)
    include_directories(${EIGEN_INCLUDE_DIRS})
    message("EIGEN_INCLUDE_DIRS" ${EIGEN_INCLUDE_DIRS})
  endif()
endif()

#add google test framework needs pthreads
find_package(GTest)
if(GTEST_FOUND)
    include_directories(${GTEST_INCLUDE_DIRS})
    set(LINK_LIBS ${LINK_LIBS} ${GTEST_BOTH_LIBRARIES})
    message("GTEST_INCLUDE_DIRS" ${GTEST_INCLUDE_DIRS})
endif()

find_package(Threads)
set(LINK_LIBS ${LINK_LIBS} ${CMAKE_THREAD_LIBS_INIT})

#define "gtest.h" include directory
#/usr/local/include/gtest/gtest.h

find_package(Doxygen)
if(DOXYGEN_FOUND)
  # Copy Doxyfile into build folder
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

  # Add custom make target 'doc'
  # Execute 'make doc' to generate documentation
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)


find_package(OpenMP)
if (OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


find_package(HDF5 COMPONENTS CXX)
if (HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIR})
  set(LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
endif()


find_package(YamlCpp)
if (YAMLCPP_FOUND)
  include_directories(${YAMLCPP_INCLUDE_DIR})
  set(LINK_LIBS "${LINK_LIBS};${YAMLCPP_LIBRARY}")
endif()


add_subdirectory(test)
add_subdirectory(examples)
add_subdirectory(data)
