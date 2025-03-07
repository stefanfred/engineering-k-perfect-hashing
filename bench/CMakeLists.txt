
add_library(kphf-bench INTERFACE)
target_link_libraries(kphf-bench PUBLIC kphf)
target_include_directory(kphf-bench PUBLIC .)
