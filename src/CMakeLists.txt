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

# find_package(SymEngine REQUIRED CONFIG) 
find_package(SymEngine REQUIRED CONFIG) 

if (SymEngine_FOUND)
    target_compile_definitions(pauliengine-objs PUBLIC HAVE_SYMENGINE)
    endif()
    
  target_link_libraries(pauliengine-objs PUBLIC symengine gmp)




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
else()
  # if you want to compile some local C++ executable for testing
  add_executable(example.x example.cpp)
  target_link_libraries(example.x PUBLIC symengine gmp)
endif()

