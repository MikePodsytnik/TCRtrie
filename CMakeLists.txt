cmake_minimum_required(VERSION 3.16)
project(tcrtrie_ LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ffast-math")

include(FetchContent)
FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.10.4
)
FetchContent_MakeAvailable(pybind11)

add_library(tcrtrie_core STATIC
        src/Trie.cpp
        src/AirrParser.cpp
)
target_include_directories(tcrtrie_core PUBLIC ${PROJECT_SOURCE_DIR}/include)
set_target_properties(tcrtrie_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

pybind11_add_module(_tcrtrie tcrtrie_/bindings.cpp)
target_link_libraries(_tcrtrie PRIVATE tcrtrie_core)

install(TARGETS _tcrtrie
        LIBRARY DESTINATION tcrtrie_
        ARCHIVE DESTINATION tcrtrie_
        RUNTIME DESTINATION tcrtrie_
)
