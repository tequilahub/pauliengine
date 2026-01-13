add_library(Math::LA INTERFACE IMPORTED)

set(BLAS_FOUND FALSE)

set(
  _blas_backends
  MKL
  BLIS
  OpenBLAS
)

foreach(backend IN LISTS _blas_backends)
  message(STATUS "Considering BLAS backend: ${backend}")

  string(TOLOWER ${backend} _backend)
  include(${CMAKE_CURRENT_LIST_DIR}/${_backend}.cmake)
  if(BLAS_FOUND)
    message(STATUS "BLAS backend ${backend} found.")
    break()
  else()
    message(STATUS "BLAS backend ${backend} not found. Trying next...")
  endif()
endforeach()

if(NOT BLAS_FOUND)
  message(FATAL_ERROR "A BLAS backend is required.")
endif()

include(CMakePrintHelpers)

cmake_print_properties(
  TARGETS
    Math::LA
  PROPERTIES
    INTERFACE_LINK_LIBRARIES
    INTERFACE_COMPILE_DEFINITIONS
    INTERFACE_COMPILE_OPTIONS
    INTERFACE_INCLUDE_DIRECTORIES
)
