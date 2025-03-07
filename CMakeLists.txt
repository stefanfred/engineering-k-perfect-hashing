
cmake_minimum_required(VERSION 3.24)
project(kPerfectHashing)

add_subdirectory(extlib/simple-ribbon)

add_library(kphf INTERFACE)
target_include_directories(kphf INTERFACE include)
target_link_libraries(kphf INTERFACE SimpleRibbon)
target_compile_features(kphf INTERFACE cxx_std_23)

add_executable(bench bench/main.cpp)
target_link_libraries(bench kphf)
target_compile_features(bench PRIVATE cxx_std_23)
