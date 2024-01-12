cmake_minimum_required(VERSION 3.27)

project(qdags)

#add_library(/Users/asugomez/lib liblsdsl libldivsufsort libldivsufsort64)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#find_library(SDSL /Users/asugomez/lib/liblsdsl)
set(GCC_COVERAGE_COMPILE_FLAGS "-std=c++11 -O3 -DNDEBUG -I /opt/homebrew/include -L /opt/homebrew/lib" )
set(CMAKE_CXX_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS}")


include_directories(/opt/homebrew/include)

separate_arguments(out UNIX_COMMAND PROGRAM SEPARATE_ARGS "~/Desktop/Magister/Tesis/qdags main.cpp ./data/prop-direct-P197 ./data/prop-direct-P800 ./data/prop-direct-P1366")

LINK_LIBRARIES(sdsl divsufsort divsufsort64)
add_executable(test2 main.cpp)

