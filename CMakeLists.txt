cmake_minimum_required(VERSION 3.24...4.0)
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
add_subdirectory(extlib/fips)

if(NOT TARGET Sux)
    add_library(Sux INTERFACE)
    target_include_directories(Sux SYSTEM INTERFACE extlib/sux)
endif()

add_library(KphfCommon INTERFACE)
target_include_directories(KphfCommon INTERFACE include/common)
target_include_directories(KphfCommon INTERFACE extlib/gcem/include)
target_compile_features(KphfCommon INTERFACE cxx_std_23)

if(NOT TARGET tlx)
    set(TLX_INSTALL_INCLUDE_DIR tlx CACHE PATH "Workaround for TLX breaking the first cmake call")
    add_subdirectory(extlib/tlx SYSTEM EXCLUDE_FROM_ALL)
endif()

#################### K-Perfect Hash Function Libraries ####################
add_library(KphfPaCHash INTERFACE)
target_include_directories(KphfPaCHash INTERFACE include/pachash)
target_link_libraries(KphfPaCHash INTERFACE ByteHamsterUtil SimpleRibbon Sux KphfCommon)
target_compile_features(KphfPaCHash INTERFACE cxx_std_23)
add_library(Kphf::PaCHash ALIAS KphfPaCHash)

add_library(KphfRecSplit INTERFACE)
target_include_directories(KphfRecSplit INTERFACE include/recsplit)
target_link_libraries(KphfRecSplit INTERFACE ByteHamsterUtil Sux KphfCommon)
target_compile_features(KphfRecSplit INTERFACE cxx_std_23)
add_library(Kphf::RecSplit ALIAS KphfRecSplit)

add_library(KphfHashDisplace INTERFACE)
target_include_directories(KphfHashDisplace INTERFACE include/hash-displace)
target_link_libraries(KphfHashDisplace INTERFACE ByteHamsterUtil KphfCommon)
target_compile_features(KphfHashDisplace INTERFACE cxx_std_23)
add_library(Kphf::HashDisplace ALIAS KphfHashDisplace)

add_library(KphfThresholdBasedBumping INTERFACE)
target_include_directories(KphfThresholdBasedBumping INTERFACE include/threshold-based-bumping)
target_link_libraries(KphfThresholdBasedBumping INTERFACE ByteHamsterUtil SimpleRibbon Sux KphfCommon FiPS::fips)
target_compile_features(KphfThresholdBasedBumping INTERFACE cxx_std_23)
add_library(Kphf::ThresholdBasedBumping ALIAS KphfThresholdBasedBumping)

#################### Benchmarks ####################
if(PROJECT_IS_TOP_LEVEL)
    add_library(PerfectHashingCompetitors INTERFACE)
    target_include_directories(PerfectHashingCompetitors INTERFACE contenders)
    target_include_directories(PerfectHashingCompetitors INTERFACE benchmarklib)
    target_link_libraries(PerfectHashingCompetitors INTERFACE ByteHamster::Util)
    target_link_libraries(PerfectHashingCompetitors INTERFACE Kphf::PaCHash Kphf::RecSplit Kphf::HashDisplace Kphf::ThresholdBasedBumping)

    file(GLOB cmphSources extlib/cmph/src/**.c)
    list(FILTER cmphSources EXCLUDE REGEX ".*cmph/src/(bdz_gen_lookup_table|bm_numbers|sdbm_hash|fnv_hash|djb2_hash|hashtree|buffer_manage|main).c$")
    add_library(Cmph SHARED ${cmphSources})
    target_include_directories(Cmph SYSTEM INTERFACE extlib/cmph/src)
    target_link_libraries(PerfectHashingCompetitors INTERFACE Cmph)
    
    add_library(PerfectHashingCompetitorParetoRunners INTERFACE)
    file(GLOB_RECURSE contendersCpp "${CMAKE_SOURCE_DIR}/contenders/**.cpp")
    foreach(filename IN LISTS contendersCpp)
        get_filename_component(name "${filename}" NAME_WLE)
        if (name STREQUAL "GpuRecSplitContender")
            continue()
        endif()
        if (name STREQUAL "GpuPhobicContender" AND NOT TARGET GpuPTHash)
            continue()
        endif()
        add_library(${name} SHARED ${filename})
        target_link_libraries(${name} PRIVATE PerfectHashingCompetitors)
        target_link_libraries(PerfectHashingCompetitorParetoRunners INTERFACE ${name})
    endforeach()

    add_executable(Comparison src/comparison.cpp)
    target_link_libraries(Comparison PRIVATE tlx ByteHamster::Util PerfectHashingCompetitors PerfectHashingCompetitorParetoRunners)
    target_compile_features(Comparison PRIVATE cxx_std_23)
endif()
