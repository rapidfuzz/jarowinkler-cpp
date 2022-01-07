cmake_minimum_required(VERSION 3.8)

project(jaro_winkler CXX)
include(GNUInstallDirs)

option(BUILD_TESTS "Build unit tests" OFF)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)

set(HEADERS
	jaro_winkler/jaro_winkler.hpp
	jaro_winkler/details/common.hpp
	jaro_winkler/details/intrinsics.hpp
	jaro_winkler/details/jaro_impl.hpp
)

add_library(jaro_winkler INTERFACE)

# provide a namespaced alias for clients to 'link' against if JaroWinkler is included as a sub-project
add_library(jaro_winkler::jaro_winkler ALIAS jaro_winkler)

set_target_properties(jaro_winkler PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_compile_features(jaro_winkler INTERFACE cxx_std_11)

target_include_directories(jaro_winkler INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
                                               "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>")

install(TARGETS jaro_winkler
  EXPORT  jaro_winkler
)

install(DIRECTORY jaro_winkler
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

if(BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif()
