cmake_minimum_required(VERSION 3.10)
project(TrieProject)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(trie_exec
    src/main.cpp
    src/Trie.cpp
)

find_package(Threads REQUIRED)
target_link_libraries(trie_exec PRIVATE Threads::Threads)
