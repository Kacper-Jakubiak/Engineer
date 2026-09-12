cmake_minimum_required(VERSION 3.21)
project(Engineer)

set(CMAKE_CXX_STANDARD 20)

add_executable(Engineer main.cpp
        include/Config.h
        include/SimulationSetup.h

        include/SchemerBuilder.h
        include/SchemerRunner.h
        include/FractionalScheme.h

        src/SchemerBuilder_setters.cpp
        src/SchemerBuilder_build.cpp
        src/SchemerRunner.cpp
        src/FractionalScheme.cpp

        capala/alfaStabilny.h
        capala/hist.h
        capala/alfaStabilny.cpp
        capala/hist.cpp
        include/Schemer.h
)

include(FetchContent)

FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 5.0.0
)

FetchContent_MakeAvailable(eigen)

target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

add_compile_definitions(-DPROJECT_ROOT="${CMAKE_SOURCE_DIR}")