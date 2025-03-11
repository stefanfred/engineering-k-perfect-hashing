cmake_minimum_required(VERSION 3.24...3.31)
project(kPerfectHashing)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release" AND PROJECT_IS_TOP_LEVEL)
    add_compile_options(-march=native)
endif()

#################### Dependency Setup ####################
add_subdirectory(extlib/simple-ribbon)

add_subdirectory(extlib/util)

if(NOT TARGET Sux)
    add_library(Sux INTERFACE)
    target_include_directories(Sux SYSTEM INTERFACE extlib/sux)
endif()

#################### K-Perfect Hash Function Libraries ####################
add_library(KphfPaCHash INTERFACE)
target_include_directories(KphfPaCHash INTERFACE include)
target_link_libraries(KphfPaCHash INTERFACE ByteHamsterUtil SimpleRibbon Sux)
target_compile_features(KphfPaCHash INTERFACE cxx_std_23)
add_library(Kphf::PaCHash ALIAS KphfPaCHash)

add_library(KphfRecSplit INTERFACE)
target_include_directories(KphfRecSplit INTERFACE include)
target_link_libraries(KphfRecSplit INTERFACE ByteHamsterUtil Sux)
target_compile_features(KphfRecSplit INTERFACE cxx_std_23)
add_library(Kphf::RecSplit ALIAS KphfRecSplit)

add_library(KphfHashDisplace INTERFACE)
target_include_directories(KphfHashDisplace INTERFACE include)
target_link_libraries(KphfHashDisplace INTERFACE ByteHamsterUtil)
target_compile_features(KphfHashDisplace INTERFACE cxx_std_23)
add_library(Kphf::HashDisplace ALIAS KphfHashDisplace)

add_library(KphfThresholdBasedBumping INTERFACE)
target_include_directories(KphfThresholdBasedBumping INTERFACE include)
target_link_libraries(KphfThresholdBasedBumping INTERFACE ByteHamsterUtil SimpleRibbon)
target_compile_features(KphfThresholdBasedBumping INTERFACE cxx_std_23)
add_library(Kphf::ThresholdBasedBumping ALIAS KphfThresholdBasedBumping)

#################### Benchmarks ####################
if(PROJECT_IS_TOP_LEVEL)
    add_executable(Benchmark bench/main.cpp)
    target_link_libraries(Benchmark ByteHamsterUtil)
    target_link_libraries(Benchmark Kphf::PaCHash Kphf::RecSplit Kphf::HashDisplace Kphf::ThresholdBasedBumping)
    target_compile_features(Benchmark PRIVATE cxx_std_23)
endif()
