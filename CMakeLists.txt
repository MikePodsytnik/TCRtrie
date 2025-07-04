cmake_minimum_required(VERSION 3.16)
project(tcrtrie LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox /fp:fast")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math")
endif()


include(FetchContent)
FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.10.4
)
FetchContent_MakeAvailable(pybind11)

add_library(tcrtrie_core STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Trie.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/src/AirrParser.cpp
)
target_include_directories(tcrtrie_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(tcrtrie_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

find_package(Threads REQUIRED)
target_link_libraries(tcrtrie_core PUBLIC Threads::Threads)

pybind11_add_module(_tcrtrie ${CMAKE_CURRENT_SOURCE_DIR}/tcrtrie/bindings.cpp)
target_link_libraries(_tcrtrie PRIVATE tcrtrie_core)

install(TARGETS _tcrtrie
        LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tcrtrie
        ARCHIVE DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tcrtrie
        RUNTIME DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/tcrtrie
)
