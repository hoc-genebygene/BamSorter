cmake_minimum_required(VERSION 3.9)
project(BamSorter CXX)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.65.1 REQUIRED filesystem)
include_directories(${Boost_INCLUDE_DIRS})

find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIRS})

find_package(LZ4 REQUIRED)
include_directories(${LZ4_INCLUDE_DIRS})

message(STATUS ${LZ4_LIBRARIES})

add_executable(BamSorter BamSorter.cpp)
target_link_libraries(BamSorter ${Boost_LIBRARIES} ${ZLIB_LIBRARIES} ${LZ4_LIBRARIES})
