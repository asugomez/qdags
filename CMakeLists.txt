cmake_minimum_required(VERSION 3.5)
project(qdags) # Defines the project name according to what we provided during project creation

set(CMAKE_CXX_STANDARD 17) #  see if change to 14 or 17

# Compiler options
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNDEBUG -DSAN=1 -fsanitize=undefined") # run with debug (assert)
#set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fsanitize=address")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG") # run with debug (assert)
include_directories(/opt/homebrew/include)
link_directories(/opt/homebrew/lib)

# main for testing
add_executable(qdags main.cpp)
# Link libraries
target_link_libraries(qdags sdsl divsufsort divsufsort64)

# ORIGINAL
add_executable(j3-original ./queries/original/j3.cpp)

# GRADUAL RETRIEVAL
add_executable(j3-partial-louds-back ./queries/partial/louds/backtracking/j3.cpp)
add_executable(j3-partial-louds-optOrder ./queries/partial/louds/optimalOrder/j3.cpp)
add_executable(j3-partial-dfuds-back ./queries/partial/dfuds/backtracking/j3.cpp)
add_executable(j3-partial-dfuds-optOrder ./queries/partial/dfuds/optimalOrder/j3.cpp)

# RANKED ENUMERATION
add_executable(j3-ranked-louds-backtracking ./queries/ranked/louds/backtracking/j3.cpp)
add_executable(j3-ranked-louds-optOrder ./queries/ranked/louds/optimalOrder/j3.cpp)
add_executable(j3-ranked-dfuds-backtracking ./queries/ranked/dfuds/backtracking/j3.cpp)
add_executable(j3-ranked-dfuds-optOrder ./queries/ranked/dfuds/optimalOrder/j3.cpp)

# LAZY QDAGS
add_executable(j3-lqdags ./queries/lqdags/j3.cpp)

# Main for testing
# ORIGINAL
add_executable(main ./main.cpp)

# Link libraries
target_link_libraries(j3-original sdsl divsufsort divsufsort64)

target_link_libraries(j3-partial-louds-back sdsl divsufsort divsufsort64)
target_link_libraries(j3-partial-louds-optOrder sdsl divsufsort divsufsort64)
target_link_libraries(j3-ranked-louds-backtracking sdsl divsufsort divsufsort64)
target_link_libraries(j3-ranked-louds-optOrder sdsl divsufsort divsufsort64)
target_link_libraries(j3-partial-dfuds-back sdsl divsufsort divsufsort64)
target_link_libraries(j3-partial-dfuds-optOrder sdsl divsufsort divsufsort64)
target_link_libraries(j3-ranked-dfuds-backtracking sdsl divsufsort divsufsort64)
target_link_libraries(j3-ranked-dfuds-optOrder sdsl divsufsort divsufsort64)

target_link_libraries(j3-lqdags sdsl divsufsort divsufsort64)

target_link_libraries(main sdsl divsufsort divsufsort64)

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


