cmake_minimum_required(VERSION 3.5)
project(qdags) # Defines the project name according to what we provided during project creation

set(CMAKE_CXX_STANDARD 14) #  see if change to 14 or 17

# Compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
include_directories(/opt/homebrew/include)
link_directories(/opt/homebrew/lib)

add_executable(qdags main.cpp)

# Link libraries
target_link_libraries(qdags sdsl divsufsort divsufsort64)

# Add command-line arguments
#add_custom_target(
#        run
#        COMMAND qdags ./data/all/prop-direct-P4673 ./data/all/prop-direct-P4144 ./data/all/prop-direct-P217 ./data/all/prop-direct-P276 0 100
#        WORKING_DIRECTORY /Users/asugomez/Desktop/Magister/Tesis/qdags
##        ./data/priorities/t4/pri1-3
##        ./data/priorities/t4/pri2-3
##        ./data/priorities/t4/pri3-3
##        ./data/priorities/t4/pri4-3
##        0
##        10000
#)


