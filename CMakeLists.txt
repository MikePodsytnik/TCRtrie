cmake_minimum_required(VERSION 3.16)
project(tcrtrie LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 CONFIG REQUIRED)

add_library(tcrtrie_core STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Trie.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AirrParser.cpp
)

target_compile_options(tcrtrie_core PRIVATE
    $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Release>>:/Ox /fp:fast>
    $<$<AND:$<NOT:$<CXX_COMPILER_ID:MSVC>>,$<CONFIG:Release>>:-O3 -ffast-math>
)

target_include_directories(tcrtrie_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(tcrtrie_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

find_package(Threads REQUIRED)
target_link_libraries(tcrtrie_core PUBLIC Threads::Threads)

pybind11_add_module(_tcrtrie ${CMAKE_CURRENT_SOURCE_DIR}/tcrtrie/bindings.cpp)
target_link_libraries(_tcrtrie PRIVATE tcrtrie_core)

install(TARGETS _tcrtrie
    LIBRARY DESTINATION tcrtrie
    ARCHIVE DESTINATION tcrtrie
    RUNTIME DESTINATION tcrtrie
)
