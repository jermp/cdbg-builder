cmake_minimum_required(VERSION 3.13)
project(cdgb-builder C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(UNIX AND (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64"))
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mno-avx512f")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2 -mpopcnt")
endif()

if(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ggdb")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-missing-braces -Wno-unknown-attributes -Wno-unused-function")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

# bits library (header-only). It bundles essentials as a sub-submodule.
set(BITS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external/bits/include)
set(ESSENTIALS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external/bits/external/essentials/include)
set(CMD_LINE_PARSER_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external/cmd_line_parser/include)
set(LZ4_DIR ${CMAKE_SOURCE_DIR}/external/lz4/lib)
set(LIBDEFLATE_DIR ${CMAKE_SOURCE_DIR}/external/libdeflate)

if(NOT EXISTS ${BITS_INCLUDE_DIR})
    message(FATAL_ERROR
        "external/bits is missing.\n"
        "Run: git submodule update --init --recursive")
endif()
if(NOT EXISTS ${CMD_LINE_PARSER_INCLUDE_DIR})
    message(FATAL_ERROR
        "external/cmd_line_parser is missing.\n"
        "Run: git submodule update --init --recursive")
endif()
if(NOT EXISTS ${LZ4_DIR}/lz4frame.h)
    message(FATAL_ERROR
        "external/lz4 is missing.\n"
        "Run: git submodule update --init --recursive")
endif()
if(NOT EXISTS ${LIBDEFLATE_DIR}/libdeflate.h)
    message(FATAL_ERROR
        "external/libdeflate is missing.\n"
        "Run: git submodule update --init --recursive")
endif()

include_directories(${BITS_INCLUDE_DIR})
include_directories(${ESSENTIALS_INCLUDE_DIR})
include_directories(${CMD_LINE_PARSER_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/external)
include_directories(${LZ4_DIR})

# LZ4 (used for intermediate bucket file compression). We only need the
# block API (LZ4_compress_default / LZ4_decompress_safe), so the static
# library is just lz4.c -- no frame, no HC, no xxhash. Custom per-spill
# framing is handled in bucket_io.hpp.
add_library(lz4_static STATIC
    ${LZ4_DIR}/lz4.c
)
set_target_properties(lz4_static PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(lz4_static PUBLIC ${LZ4_DIR})

# libdeflate (decompress-only). Vendored as a submodule and built via
# its own CMakeLists.txt with compression / zlib / shared lib / cli
# disabled. We only call libdeflate_alloc_decompressor +
# libdeflate_gzip_decompress_ex from seq_reader.hpp -- no compressor
# state, no zlib wrapper. ~2-3x faster than zlib for decompression on
# bacterial genome FASTA inputs.
set(LIBDEFLATE_BUILD_STATIC_LIB ON CACHE BOOL "" FORCE)
set(LIBDEFLATE_BUILD_SHARED_LIB OFF CACHE BOOL "" FORCE)
set(LIBDEFLATE_COMPRESSION_SUPPORT OFF CACHE BOOL "" FORCE)
set(LIBDEFLATE_DECOMPRESSION_SUPPORT ON CACHE BOOL "" FORCE)
set(LIBDEFLATE_GZIP_SUPPORT ON CACHE BOOL "" FORCE)
set(LIBDEFLATE_ZLIB_SUPPORT OFF CACHE BOOL "" FORCE)
set(LIBDEFLATE_BUILD_GZIP OFF CACHE BOOL "" FORCE)
set(LIBDEFLATE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LIBDEFLATE_INSTALL OFF CACHE BOOL "" FORCE)
add_subdirectory(${LIBDEFLATE_DIR} ${CMAKE_BINARY_DIR}/libdeflate EXCLUDE_FROM_ALL)

find_package(Threads REQUIRED)

add_executable(cdgb-build src/main.cpp)
target_link_libraries(cdgb-build PRIVATE Threads::Threads lz4_static libdeflate_static)
target_include_directories(cdgb-build PRIVATE ${LIBDEFLATE_DIR})

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS}")
