cmake_minimum_required(VERSION 3.21)
project(Engineer)

set(CMAKE_CXX_STANDARD 20)

add_executable(Engineer main.cpp)

include(FetchContent)

FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 5.0.0
)

FetchContent_MakeAvailable(eigen)

target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)