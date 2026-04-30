cmake_minimum_required(VERSION 3.13)
project(cdgb-builder CXX)

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

if(NOT EXISTS ${BITS_INCLUDE_DIR})
    message(FATAL_ERROR
        "external/bits is missing.\n"
        "Run: git submodule update --init --recursive")
endif()

include_directories(${BITS_INCLUDE_DIR})
include_directories(${ESSENTIALS_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/external)

# zlib for transparent .gz input via kseq.h
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

add_executable(cdgb-build src/main.cpp)
target_link_libraries(cdgb-build PRIVATE ZLIB::ZLIB Threads::Threads)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ flags: ${CMAKE_CXX_FLAGS}")
