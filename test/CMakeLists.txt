# Test executable. Include dirs (bits, essentials, unordered_dense,
# include/, external/, lz4) are inherited from the parent CMakeLists via
# include_directories(). It links the same static libs the main binary does.

add_executable(test_stitch test_stitch.cpp)
target_link_libraries(test_stitch PRIVATE Threads::Threads lz4_static libdeflate_static)
target_include_directories(test_stitch PRIVATE ${LIBDEFLATE_DIR})
