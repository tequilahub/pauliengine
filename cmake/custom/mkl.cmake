set(MKL_INTERFACE "lp64")
set(MKL_LINK "dynamic")
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(MKL_THREADING "gnu_thread")
else()
  set(MKL_THREADING "intel_thread")
endif()

# attempt to find MKL using the provided CMake module
find_package(MKL QUIET)
if(TARGET MKL::MKL)
  message(STATUS "Using MKL as linear algebra backend")

  target_compile_definitions(
    Math::LA
    INTERFACE
      "USE_MKL"
      "MKL_LP64"
  )
  target_compile_options(
    Math::LA
    INTERFACE
      $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>
  )
  target_include_directories(
    Math::LA
    INTERFACE
      $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>
  )
  target_link_libraries(Math::LA INTERFACE $<LINK_ONLY:MKL::MKL>)

  set(BLAS_FOUND TRUE)
endif()
