# attempt to find OpenBLAS using the provided CMake module
find_package(
  OpenBLAS
  QUIET
  COMPONENTS
    shared
    openmp
)

if(TARGET OpenBLAS::OpenBLAS)
  # if the module defines a target, use it
  message(STATUS "Using OpenBLAS as linear algebra backend")

  # remove OpenMP components from Math::LA interface link libraries, if present.
  # this is a hack to avoid propagating OpenMP_C and OpenMP_Fortran from the OpenBLAS target
  get_target_property(
    _math_la_interface_libs
    OpenBLAS::OpenBLAS
    INTERFACE_LINK_LIBRARIES
  )

  if(_math_la_interface_libs)
    set(_filtered_math_la_interface_libs "")
    foreach(_lib IN LISTS _math_la_interface_libs)
      if(
        _lib
          STREQUAL
          "OpenMP::OpenMP_C"
        OR
          _lib
            STREQUAL
            "OpenMP::OpenMP_Fortran"
      )
        continue()
      endif()
      list(APPEND _filtered_math_la_interface_libs "${_lib}")
    endforeach()

    if(NOT _filtered_math_la_interface_libs STREQUAL _math_la_interface_libs)
      set_property(
        TARGET
          OpenBLAS::OpenBLAS
        PROPERTY
          INTERFACE_LINK_LIBRARIES
            "${_filtered_math_la_interface_libs}"
      )
      message(
        STATUS
        "Stripped OpenMP::OpenMP_C and OpenMP::OpenMP_Fortran from OpenBLAS::OpenBLAS"
      )
    endif()
  endif()

  target_compile_definitions(Math::LA INTERFACE "USE_OPENBLAS")
  target_compile_options(
    Math::LA
    INTERFACE
      $<TARGET_PROPERTY:OpenBLAS::OpenBLAS,INTERFACE_COMPILE_OPTIONS>
  )
  target_include_directories(
    Math::LA
    INTERFACE
      $<TARGET_PROPERTY:OpenBLAS::OpenBLAS,INTERFACE_INCLUDE_DIRECTORIES>
  )
  target_link_libraries(Math::LA INTERFACE $<LINK_ONLY:OpenBLAS::OpenBLAS>)

  set(BLAS_FOUND TRUE)
elseif(OpenBLAS_LIBRARIES AND OpenBLAS_INCLUDE_DIRS)
  # the module does not define a target, but provides variables
  message(
    STATUS
    "Using OpenBLAS as linear algebra backend (found via variables)"
  )

  target_compile_definitions(Math::LA INTERFACE "USE_OPENBLAS")
  target_include_directories(Math::LA INTERFACE ${OpenBLAS_INCLUDE_DIRS})
  target_link_libraries(Math::LA INTERFACE ${OpenBLAS_LIBRARIES})

  set(BLAS_FOUND TRUE)
else()
  # the module likely does not exist, use FindBLAS
  set(BLA_VENDOR "OpenBLAS")
  find_package(BLAS QUIET)
  message(
    STATUS
    "Using OpenBLAS as linear algebra backend (found via FindBLAS.cmake)"
  )

  # Extract the include directory from BLAS_LIBRARIES path
  get_filename_component(_blas_lib_dir "${BLAS_LIBRARIES}" DIRECTORY)
  set(_blas_include_dir "${_blas_lib_dir}")

  if(_blas_include_dir MATCHES "(/|^)(lib64)(/|$)")
    string(
      REGEX REPLACE
      "(/|^)lib64(/|$)"
      "\\1include\\2"
      _blas_include_dir
      "${_blas_include_dir}"
    )
  elseif(_blas_include_dir MATCHES "(/|^)(lib)(/|$)")
    string(
      REGEX REPLACE
      "(/|^)lib(/|$)"
      "\\1include\\2"
      _blas_include_dir
      "${_blas_include_dir}"
    )
  endif()

  include(${CMAKE_CURRENT_LIST_DIR}/cblas-lapacke.cmake)

  set(CMAKE_REQUIRED_QUIET OFF)

  # get directory containing cblas.h
  # we are assuming this is on an already-known system path
  find_path(
    _include_dirs
    NAMES
      cblas.h
    HINTS
      ${CMAKE_PREFIX_PATH}/include
    PATHS
      ${_blas_include_dir}/openblas
      ${_blas_include_dir}/openblas-openmp
      ${_blas_include_dir}
  )
  list(APPEND _la_include_dirs ${_include_dirs})

  set(CMAKE_REQUIRED_QUIET ON)

  target_compile_definitions(Math::LA INTERFACE "USE_OPENBLAS")
  target_include_directories(Math::LA INTERFACE ${_la_include_dirs})
  target_link_libraries(Math::LA INTERFACE ${BLAS_LIBRARIES})

  set(BLAS_FOUND TRUE)
endif()
