# transform pauliengine_CXX_FLAGS and EXTRA_CXXFLAGS to ;-separated lists
message(STATUS "pauliengine_CXX_FLAGS : ${pauliengine_CXX_FLAGS}")
if(DEFINED pauliengine_CXX_FLAGS)
  string(REPLACE " " ";" pauliengine_CXX_FLAGS ${pauliengine_CXX_FLAGS})
endif()
if(DEFINED EXTRA_CXXFLAGS)
  string(REPLACE " " ";" EXTRA_CXXFLAGS ${EXTRA_CXXFLAGS})
endif()

# create object library for core pauliengine functionality
add_library(
  pauliengine-objs 
  OBJECT
  QubitHamiltonian.cpp)

# target_include_directories(pauliengine-objs INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(
  pauliengine-objs
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

# No external dependencies: symbolic coefficients are the header-only polynomials
# of include/pauliengine/Polynomial.h.

if(OpenMP_CXX_FOUND)
    target_link_libraries(pauliengine-objs PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(pauliengine-objs PUBLIC PAULIENGINE_HAS_OPENMP)
endif()
    
# flags to prepend
target_compile_options(
  pauliengine-objs
  BEFORE
  PUBLIC
    "${pauliengine_CXX_FLAGS}"
    "${ARCH_FLAG}"
)

# flags to append
target_compile_options(pauliengine-objs PUBLIC "${EXTRA_CXXFLAGS}")

if(DEFINED SKBUILD)
  add_subdirectory(bindings)
endif()

