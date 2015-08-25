cmake_minimum_required(VERSION 2.8)

project(waveblocks)

add_definitions(
    -std=c++11
    -Wall
    -Wno-deprecated-declarations
    -Wextra
    -pedantic
)

add_library(waveblocks-static STATIC 
    waveblocks/csv/coefficients_file_parser.cpp
    waveblocks/util/complexnumber_parser.cpp
)

# needs pkg-config
include(FindPkgConfig)
if(PKG_CONFIG_FOUND)
    # pkg_search_module(<XPREFIX> [REQUIRED] [QUIET] <MODULE> [<MODULE>]*)
    # checks for given modules and uses the first working one
    #
    # NOTE: <XPREFIX> is case sensitive!
    #
    #  @return: PKG_CONFIG_FOUND ... if pkg-config executable was found
    #  @return: PKG_CONFIG_EXECUTABLE ... pathname of the pkg-config program
    #  @return: PKG_CONFIG_VERSION_STRING ... the version of the pkg-config program found
    #  @return: <XPREFIX>_FOUND ... set to 1 if module(s) exist
    #  @return: <XPREFIX>_LIBRARIES ... only the libraries (w/o the '-l')
    #  @return: <XPREFIX>_LIBRARY_DIRS ... the paths of the libraries (w/o the '-L')
    #  @return: <XPREFIX>_LDFLAGS ... all required linker flags
    #  @return: <XPREFIX>_LDFLAGS_OTHER ... all other linker flags
    #  @return: <XPREFIX>_INCLUDE_DIRS ... the '-I' preprocessor flags (w/o the '-I')
    #  @return: <XPREFIX>_CFLAGS ... all required cflags
    #  @return:<XPREFIX>_CFLAGS_OTHER ... the other compiler flags
    pkg_search_module(EIGEN REQUIRED eigen3)
    
    if(EIGEN_FOUND)
        # target_include_directories(waveblocks PUBLIC ${EIGEN_INCLUDE_DIRS})
        include_directories(${EIGEN_INCLUDE_DIRS})
    endif()
endif()


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

add_subdirectory(test)
add_subdirectory(examples)
add_subdirectory(waveblocks/misc)