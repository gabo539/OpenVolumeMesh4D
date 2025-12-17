cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
cmake_policy(VERSION 3.12.0)


add_executable(generate_shapes generate_shapes.cc
        subdivision.hh)

target_link_libraries(generate_shapes OpenVolumeMesh::OpenVolumeMesh)

set_target_properties(generate_shapes PROPERTIES
        CXX_STANDARD ${OVM_CXX_STANDARD}
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)