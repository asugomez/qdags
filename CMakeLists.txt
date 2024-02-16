cmake_minimum_required(VERSION 3.1...3.28)
project(
        rankedQdagsProject
        VERSION 1.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

# Añade la ruta de búsqueda de encabezados de la librería
include_directories(/opt/homebrew/include/sdsl)
include_directories(/opt/homebrew/include/divsufsort)
include_directories(/opt/homebrew/include/divsufsort64)

#add_library(sdsl STATIC IMPORTED)
#add_library(divsufsort STATIC IMPORTED)
#add_library(divsufsort64 STATIC IMPORTED)


#add_library(LibBits simple_lib.cpp simple_lib.hpp)
add_executable(main main.cpp
        src/qdagWeight.cpp
        src/joins_ranked_size.cpp)

add_library(utils src/utils.cpp)
