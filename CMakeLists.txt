cmake_minimum_required(VERSION 3.5)
project(
        rankedQdagsProject
        VERSION 1.0
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

add_executable(main main.cpp)

# Lista de archivos fuente dentro de la carpeta src
set(SOURCE_FILES
        src/joins.cpp
        src/joins_ranked_size.cpp
        src/qdags.hpp
        src/rank.hpp
        src/se_quadtree.hpp
)

# Agrega los archivos fuente al ejecutable principal
target_sources(main PRIVATE ${SOURCE_FILES})

# Añade la ruta de búsqueda de encabezados de la librería
target_link_libraries(main PRIVATE sdsl divsufsort divsufsort64)
target_include_directories(main PRIVATE /opt/homebrew/include)
target_link_directories(main PRIVATE /opt/homebrew/lib)

target_compile_options(main PRIVATE -O3 -O0 -DNDEBUG)

