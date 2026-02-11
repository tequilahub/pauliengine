set(BLA_VENDOR "FLAME")
find_package(BLAS QUIET)

if(TARGET BLAS::BLAS)
  message(
    STATUS
    "Using BLIS as linear algebra backend (found via FindBLAS.cmake)"
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
      blis.h
    HINTS
      ${CMAKE_PREFIX_PATH}/include
    PATHS
      ${_blas_include_dir}/blis
      ${_blas_include_dir}/blis-openmp
      ${_blas_include_dir}
  )
  list(APPEND _la_include_dirs ${_include_dirs})

  set(CMAKE_REQUIRED_QUIET ON)

  target_compile_definitions(Math::LA INTERFACE "USE_BLIS")
  target_include_directories(Math::LA INTERFACE ${_la_include_dirs})
  target_link_libraries(Math::LA INTERFACE ${BLAS_LIBRARIES})

  set(BLAS_FOUND TRUE)
endif()
