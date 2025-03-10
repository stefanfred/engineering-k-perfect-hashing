cmake_minimum_required(VERSION 3.24)
project(kPerfectHashing)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release" AND PROJECT_IS_TOP_LEVEL)
    add_compile_options(-march=native)
endif()

add_subdirectory(extlib/simple-ribbon)

add_subdirectory(extlib/util)

add_library(kphf INTERFACE)
target_include_directories(kphf INTERFACE include)
target_link_libraries(kphf INTERFACE SimpleRibbon ByteHamsterUtil)
target_compile_features(kphf INTERFACE cxx_std_23)

if(PROJECT_IS_TOP_LEVEL)
    add_executable(bench bench/main.cpp)
    target_link_libraries(bench kphf)
    target_compile_features(bench PRIVATE cxx_std_23)
endif()
