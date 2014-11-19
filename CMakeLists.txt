cmake_minimum_required(VERSION 2.8)

# Find the BEM++ library installed on your disk
find_package(Bempp REQUIRED CONFIG)
find_library(BEMPP_LIBRARY bempp PATHS ${BEMPP_LIBRARY_DIR})
# Add the BEM++ installation directory to CMake search path to ensure that
# e.g. the version of SWIG coming with BEM++ is used
set(CMAKE_PREFIX_PATH "${BEMPP_PREFIX}/bempp;${CMAKE_PREFIX_PATH}"
    CACHE STRING "")
# Make the new module be installed in the main BEM++ directory, by default
set(CMAKE_INSTALL_PREFIX "${BEMPP_PREFIX}/bempp" CACHE PATH "" FORCE)
# By default, compile the new module with optimizations turned on
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
       FORCE)
endif ()

# The simple_vector_spaces library
add_library(simple_vector_spaces SHARED 
    simple_vector_space.cpp
    piecewise_constant_vector_space.cpp 
    piecewise_linear_vector_space.cpp
    piecewise_linear_continuous_vector_space.cpp
    piecewise_linear_discontinuous_vector_space.cpp
)
target_link_libraries(simple_vector_spaces ${BEMPP_LIBRARY} ${BEMPP_TEUCHOS_LIBRARY})
set_target_properties(simple_vector_spaces PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/bempp/lib")
install(TARGETS simple_vector_spaces LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/bempp/lib")

# The integrate_grid_function library
add_library(integrate_grid_function SHARED 
    integrate_grid_function.cpp
)
target_link_libraries(integrate_grid_function ${BEMPP_LIBRARY} ${BEMPP_TEUCHOS_LIBRARY})
set_target_properties(integrate_grid_function PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/bempp/lib")
install(TARGETS integrate_grid_function LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/bempp/lib")

# Find SWIG

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})
set(CMAKE_SWIG_FLAGS "-modern")

# Add include directories necessary for compilation of the new modules

include_directories(${BEMPP_PYTHON_INCLUDE_DIRS}) # Python and NumPy headers
include_directories("${BEMPP_INCLUDE_DIR}")
include_directories("${BEMPP_INCLUDE_DIR}/bempp")
include_directories("${BEMPP_INCLUDE_DIR}/bempp/swig")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

# Specify the source files and name of the new modules

set_source_files_properties(simple_vector_spaces.i PROPERTIES CPLUSPLUS ON)
swig_add_module(simple_vector_spaces python simple_vector_spaces.i)
swig_link_libraries(simple_vector_spaces
    simple_vector_spaces
    ${BEMPP_PYTHON_LIBRARY}
    ${BEMPP_LIBRARY}
    ${BEMPP_TEUCHOS_LIBRARY})

set_source_files_properties(integrate_grid_function.i PROPERTIES CPLUSPLUS ON)
swig_add_module(integrate_grid_function python integrate_grid_function.i)
swig_link_libraries(integrate_grid_function
    integrate_grid_function
    ${BEMPP_PYTHON_LIBRARY}
    ${BEMPP_LIBRARY}
    ${BEMPP_TEUCHOS_LIBRARY})

# Specify files to be installed
set(FILES_TO_INSTALL
    "${CMAKE_BINARY_DIR}/simple_vector_spaces.py"
    "${CMAKE_BINARY_DIR}/_simple_vector_spaces.so")
install(FILES ${FILES_TO_INSTALL}
    DESTINATION "${CMAKE_INSTALL_PREFIX}/python/bempp")

# Specify files to be installed
set(FILES_TO_INSTALL
    "${CMAKE_BINARY_DIR}/integrate_grid_function.py"
    "${CMAKE_BINARY_DIR}/_integrate_grid_function.so")
install(FILES ${FILES_TO_INSTALL}
    DESTINATION "${CMAKE_INSTALL_PREFIX}/python/bempp")
