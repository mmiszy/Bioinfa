cmake_minimum_required(VERSION 3.3)
project(Bioinfa)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    Alignment.cpp
    AlignmentTest.cpp
    main.cpp Alignment.h)

add_executable(Bioinfa ${SOURCE_FILES})