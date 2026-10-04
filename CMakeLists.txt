cmake_minimum_required(VERSION 3.16)
project(tcrtrie LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)

find_package(Threads REQUIRED)

add_executable(tcrtrie
        main/main.cpp
        src/CliConfig.cpp
        src/CostMatrix.cpp
        src/TsvReader.cpp
        src/Trie.cpp
        src/RegionTrie.cpp
        src/MatchWriter.cpp
        src/RepertoireMatcher.cpp
)

target_compile_options(tcrtrie PRIVATE
        $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Release>>:/Ox /fp:fast>
        $<$<AND:$<NOT:$<CXX_COMPILER_ID:MSVC>>,$<CONFIG:Release>>:-O3 -ffast-math>
)

target_include_directories(tcrtrie PRIVATE include)
target_link_libraries(tcrtrie PRIVATE Threads::Threads)

install(TARGETS tcrtrie
        RUNTIME DESTINATION tcrtrie
        COMPONENT python
)
