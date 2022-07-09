cmake_minimum_required (VERSION 3.16)
project(manif-geom-cpp)

set(CMAKE_CXX_STANDARD 17)

option(BUILD_TESTS "Build Tests" ON)

find_package(Eigen3 REQUIRED NO_MODULE)

add_library(manif-geom-cpp INTERFACE)
target_include_directories(manif-geom-cpp INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(manif-geom-cpp INTERFACE Eigen3::Eigen)

set(UNIT_TEST unit-tests)
add_executable(${UNIT_TEST}
    tests/MainTest.cpp
    tests/SO2Test.cpp
    tests/SE2Test.cpp
    tests/SO3Test.cpp
    tests/SE3Test.cpp
)
target_link_libraries(${UNIT_TEST}
    manif-geom-cpp
    boost_unit_test_framework
)
add_test(NAME ${UNIT_TEST} COMMAND ${UNIT_TEST})

add_custom_command(
    TARGET ${UNIT_TEST}
    COMMENT "Run manifold unit tests"
    POST_BUILD
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND ${UNIT_TEST}
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/manif-geom-cppConfigVersion.cmake"
    VERSION 0.0
    COMPATIBILITY AnyNewerVersion
)

install(TARGETS manif-geom-cpp
    EXPORT manif-geom-cppTargets
    LIBRARY DESTINATION lib COMPONENT Runtime
    ARCHIVE DESTINATION lib COMPONENT Development
    RUNTIME DESTINATION bin COMPONENT Runtime
    PUBLIC_HEADER DESTINATION include COMPONENT Development
    BUNDLE DESTINATION bin COMPONENT Runtime
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/manif-geom-cppConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/manif-geom-cppConfig.cmake"
    INSTALL_DESTINATION lib/cmake/manif-geom-cpp
)

install(EXPORT manif-geom-cppTargets DESTINATION lib/cmake/manif-geom-cpp)
install(FILES "${PROJECT_BINARY_DIR}/manif-geom-cppConfigVersion.cmake"
              "${PROJECT_BINARY_DIR}/manif-geom-cppConfig.cmake"
        DESTINATION lib/cmake/manif-geom-cpp)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)