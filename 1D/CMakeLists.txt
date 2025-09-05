cmake_minimum_required(VERSION 3.10)
project(chaleur_1D)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

add_executable(chaleur_1D
    src/main.cpp
    src/solver.cpp
    src/boundary.cpp
    src/source_term.cpp
)
