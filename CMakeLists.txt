cmake_minimum_required(VERSION 3.21)
project(Engineer)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
        Schemer
        GIT_REPOSITORY https://github.com/Kacper-Jakubiak/Schemer.git
        GIT_TAG        master
)
FetchContent_MakeAvailable(Schemer)

add_executable(Engineer
        main.cpp
        capala/alfaStabilny.h
        capala/hist.h
        capala/alfaStabilny.cpp
        capala/hist.cpp
)

target_link_libraries(Engineer PRIVATE Schemer::Schemer)

target_compile_definitions(Engineer PUBLIC -DPROJECT_ROOT="${CMAKE_SOURCE_DIR}")