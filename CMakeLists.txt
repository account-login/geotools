set(CMAKE_LEGACY_CYGWIN_WIN32 0)
cmake_minimum_required(VERSION 2.8)

project(geotools)

set(CMAKE_CXX_STANDARD 11)

add_definitions(-Wall -Wextra -g -fno-access-control)
add_definitions(-O2)

add_library(DoctestMain OBJECT tests/doctest_main.cpp)

add_executable(bench_query
    tests/bench_query.cpp

    # for clion
    lruset.hpp
    geoutil.hpp
    geotree.hpp
    geodensity.hpp
    geodensity_bounded.hpp
    edgetree.cpp
)
target_link_libraries(bench_query
    boost_system
)

add_executable(test_lruset
    tests/test_lruset.cpp
    $<TARGET_OBJECTS:DoctestMain>
)

add_executable(test_geoutil
    tests/test_geoutil.cpp
    $<TARGET_OBJECTS:DoctestMain>
)

add_executable(test_geotree
    tests/test_geotree.cpp
    $<TARGET_OBJECTS:DoctestMain>
)

add_executable(test_edgetree
    tests/test_edgetree.cpp
    $<TARGET_OBJECTS:DoctestMain>
)

add_executable(test_geodensity
    tests/test_geodensity.cpp
    $<TARGET_OBJECTS:DoctestMain>
)
target_link_libraries(test_geodensity
    boost_thread boost_system
)

add_executable(test_geodensity_bounded
    tests/test_geodensity_bounded.cpp
    $<TARGET_OBJECTS:DoctestMain>
)
target_link_libraries(test_geodensity_bounded
    boost_thread boost_system
)
