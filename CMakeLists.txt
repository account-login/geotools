set(CMAKE_LEGACY_CYGWIN_WIN32 0)
cmake_minimum_required(VERSION 2.8)

project(geotools)

set(CMAKE_CXX_STANDARD 11)

add_executable(bench_query
    tests/bench_query.cpp
    tests/catch.cpp
    # for clion
    lruset.hpp
    geoutil.hpp
    geotree.hpp
    geodensity.hpp
    geodensity_bounded.hpp
)

add_executable(test_lruset
    tests/test_lruset.cpp
    tests/catch.cpp
)

add_executable(test_geoutil
    tests/test_geoutil.cpp
    tests/catch.cpp
)

add_executable(test_geotree
    tests/test_geotree.cpp
    tests/catch.cpp
)

add_executable(test_geodensity
    tests/test_geodensity.cpp
    tests/catch.cpp
)
target_link_libraries(test_geodensity
    boost_thread boost_system
)

add_executable(test_geodensity_bounded
    tests/test_geodensity_bounded.cpp
    tests/catch.cpp
)
target_link_libraries(test_geodensity_bounded
    boost_thread boost_system
)
