# qdags
 
To compile:
```
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib program.cpp -o program -lsdsl -ldivsufsort -ldivsufsort64 
```
-std=c++11 -O3 -DNDEBUG -I /opt/homebrew/include -L /opt/homebrew/lib -lsdsl -ldivsufsort -ldivsufsort64

```
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64 program.cpp -o program  
```

To run main 
```
./main ./data/prop-direct-P197 ./data/prop-direct-P800 ./data/prop-direct-P1366 
```

/usr/local/include --> divsufsort ,...
See http://algo2.iti.kit.edu/gog/docs/html/index.html 

# set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message("CXX Standard: ${CMAKE_CXX_STANDARD}")

# variables to store additional compiler flags
set(CXX_FLAGS)
set(CXX_FLAGS_DEBUG)
set(CXX_FLAGS_RELEASE)

# always build shared libraries by default
set(BUILD_SHARED_LIBS ON)

#add_compile_options("-std=c++11")

# append global compiler flags to list
list(APPEND CXX_FLAGS "-std=c++11" "-O3" "-DNDEBUG" "-lsdsl" "-ldivsufsort" "-ldivsufsort64" )
CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -L$(LIB_DIR)
CCLIB=-lsdsl -ldivsufsort -ldivsufsort64
LIB_DIR = @CMAKE_INSTALL_PREFIX@/lib
INC_DIR = @CMAKE_INSTALL_PREFIX@/include
MY_CXX_FLAGS=@CMAKE_CXX_FLAGS@ $(CODE_COVER)
MY_CXX_OPT_FLAGS=@CMAKE_CXX_OPT_FLAGS@
MY_CXX=@CMAKE_CXX_COMPILER@
MY_CC=@CMAKE_C_COMPILER@

#set(GCC_COVERAGE_COMPILE_FLAGS "-O3 -DNDEBUG -I" )
#set(GCC_COVERAGE_LINK_FLAGS    "-lsdsl -ldivsufsort -ldivsufsort64")

#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
#SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")


# output compiler flags
include(CMakePrintHelpers)
cmake_print_variables(CXX_FLAGS)
cmake_print_variables(CXX_FLAGS_DEBUG)
cmake_print_variables(CXX_FLAGS_RELEASE)


# set compile options for interp library
target_compile_options(qdags PRIVATE ${CXX_FLAGS})



# set compile options for the main executable
target_compile_options(main PRIVATE ${CXX_FLAGS}
"$<$<CONFIG:Debug>:${CXX_FLAGS_DEBUG}>"
"$<$<CONFIG:Release>:${CXX_FLAGS_RELEASE}>")

# link to shared library
target_link_libraries(main qdags)

# generate configuration file
#configure_file(config.h.in config.h @ONLY)

# include binary  directory to find config.h
target_include_directories(main PRIVATE ${PROJECT_BINARY_DIR})

# compile an object library for producing a shared and dynamic library
target_link_directories(qdags /Users/asugomez/lib/)
target_link_libraries(qdags sdsl)
add_library(qdagsLib
/Users/asugomez/lib/libsdsl.a
/Users/asugomez/lib/libdivsufsort.a
/Users/asugomez/lib/libdivsufsort64.a
./src/joins.cpp
)


# Tests
In the queries folder we have all the queries to be tested grouped by the type of test: 
- partial: 
  - dfuds:
    - fixedQueue
    - nonFixedQueue
  - louds:
    - fixedQueue
    - nonFixedQueue
- ranked:
    - dfuds:
        - fixedQueue
        - nonFixedQueue
    - louds:
        - fixedQueue
        - nonFixedQueue
- lqdags

The script files to run the tests are in the runqueries folder. On each subfolder we have:
- Compile script: it compiles all the queries.
- Run script: it runs all the queries and put the outputs in the respective output folder.
- All the run queries to be tested with different data. The files are named as runqueries-[type of query]-bfs-sorted.sh

## To run the tests
1. First go to the runqueries folder and choose the subfolder you want to test. For example
` cd runqueries/partial/dfuds/fixedQueue`
2. Run the compile script
`./compilePartialDfudsFixedQueue.sh`
3. Run the run script
`./runPartialDfudsFixed.sh`