# transform PauliEngine_CXX_FLAGS and EXTRA_CXXFLAGS to ;-separated lists
message(STATUS "PauliEngine_CXX_FLAGS : ${PauliEngine_CXX_FLAGS}")
if(DEFINED PauliEngine_CXX_FLAGS)
  string(REPLACE " " ";" PauliEngine_CXX_FLAGS ${PauliEngine_CXX_FLAGS})
endif()
if(DEFINED EXTRA_CXXFLAGS)
  string(REPLACE " " ";" EXTRA_CXXFLAGS ${EXTRA_CXXFLAGS})
endif()

# create object library for core PauliEngine functionality
add_library(pauli_string INTERFACE)

target_include_directories(pauli_string INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(
#   pauli_string
#   PUBLIC
#     $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
#     $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
# )

# find_package(SymEngine REQUIRED CONFIG) 
find_package(SymEngine CONFIG) 

if (SymEngine_FOUND)
    target_compile_definitions(pauli_string INTERFACE HAVE_SYMENGINE)
    target_link_libraries(pauli_string INTERFACE symengine)
endif()



# flags to prepend
target_compile_options(
  pauli_string
  BEFORE
  INTERFACE
    "${PauliEngine_CXX_FLAGS}"
    "${ARCH_FLAG}"
)

# target_link_libraries(
#   pauli_string
#   INTERFACE
#     Math::LA
#     OpenMP::OpenMP_CXX
# )

# flags to append
target_compile_options(pauli_string INTERFACE "${EXTRA_CXXFLAGS}")

if(DEFINED SKBUILD)
  add_subdirectory(bindings)
else()
  # if you want to compile some local C++ executable for testing
  #add_executable(example.x example.cpp)
  #target_link_libraries(example.x PRIVATE pauli_string)
endif()