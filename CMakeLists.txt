# Cmake config largely taken from catch2
cmake_minimum_required(VERSION 3.5)

# detect if Catch is being bundled,
# disable testsuite in that case
if(NOT DEFINED PROJECT_NAME)
  set(NOT_SUBPROJECT ON)
else()
  set(NOT_SUBPROJECT OFF)
endif()

option(JARO_WINKLER_BUILD_TESTING "Build tests" OFF)
option(JARO_WINKLER_BUILD_BENCHMARKS "Build benchmarks" OFF)

# JaroWinkler's build breaks if done in-tree. You probably should not build
# things in tree anyway, but we can allow projects that include JaroWinkler
# as a subproject to build in-tree as long as it is not in our tree.
if (CMAKE_BINARY_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    message(FATAL_ERROR "Building in-source is not supported! Create a build dir and remove ${CMAKE_SOURCE_DIR}/CMakeCache.txt")
endif()

project(JaroWinkler LANGUAGES CXX VERSION 0.0.1)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Basic paths
set(BASE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(SOURCES_DIR ${BASE_DIR}/jaro_winkler)
set(TEST_DIR ${BASE_DIR}/test)
set(BENCHMARK_DIR ${BASE_DIR}/tests/bench)
set(EXAMPLES_DIR ${BASE_DIR}/examples)

# Please keep these ordered alphabetically
set(PUBLIC_HEADERS
    ${SOURCES_DIR}/jaro_winkler.hpp
)

add_library(jaro_winkler INTERFACE)

# provide a namespaced alias for clients to 'link' against if JaroWinkler is included as a sub-project
add_library(jaro_winkler::jaro_winkler ALIAS jaro_winkler)

set_target_properties(jaro_winkler PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")
target_compile_features(jaro_winkler INTERFACE cxx_std_14)

target_include_directories(jaro_winkler
    INTERFACE
      $<BUILD_INTERFACE:${SOURCES_DIR}/..>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Build tests only if requested
if(JARO_WINKLER_BUILD_TESTING AND NOT_SUBPROJECT)
    include(CTest)
    enable_testing()
    add_subdirectory(test)
endif()

# Build examples only if requested
if(JARO_WINKLER_BUILD_EXAMPLES)
    # todo write examples
    #add_subdirectory(examples)
endif()

# Build benchmarks only if requested
if(JAROW_WINKLER_BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif()

# Build fuzz tests only if requested
if(JARO_WINKLER_BUILD_FUZZERS)
    # todo write fuzz tests
    # add_subdirectory(fuzzing)
endif()

# Only perform the installation steps when JaroWinkler is not being used as
# a subproject via `add_subdirectory`
if (NOT_SUBPROJECT)
    set(JARO_WINKLER_CMAKE_CONFIG_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/rapidfuzz")

    install(
        TARGETS
          jaro_winkler
        EXPORT
          JaroWinklerTargets
        DESTINATION
          ${CMAKE_INSTALL_LIBDIR}
    )

    install(
        EXPORT
          JaroWinklerTargets
        NAMESPACE
          jaro_winkler::
        DESTINATION
          ${JARO_WINKLER_CMAKE_CONFIG_DESTINATION}
    )

    install(
        DIRECTORY
          jaro_winkler
        DESTINATION
          ${CMAKE_INSTALL_INCLUDEDIR}
        FILES_MATCHING
          PATTERN "*.hpp"
          PATTERN "*.impl"
    )

    configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        INSTALL_DESTINATION ${JARO_WINKLER_CMAKE_CONFIG_DESTINATION}
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        COMPATIBILITY SameMajorVersion
    )

    install(
      FILES
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
      DESTINATION
        ${JARO_WINKLER_CMAKE_CONFIG_DESTINATION}
    )

    # CPack/CMake started taking the package version from project version 3.12
    # So we need to set the version manually for older CMake versions
    if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
        set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
    endif()
    
    set(CPACK_PACKAGE_VENDOR "Max Bachmann")
    set(CPACK_PACKAGE_CONTACT "https://github.com/maxbachmann/jarowinkler-cpp")
    include(CPack)

endif(NOT_SUBPROJECT)