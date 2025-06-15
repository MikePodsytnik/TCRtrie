cmake_minimum_required(VERSION 3.10)
project(TrieProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math")

include_directories(${PROJECT_SOURCE_DIR}/include)

include(FetchContent)
FetchContent_Declare(
        CLI11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG v2.3.2
)
FetchContent_MakeAvailable(CLI11)

add_executable(TCRtrie
        src/main.cpp
        src/Trie.cpp
        src/TrieInterface.cpp
        src/AirrParser.cpp
)

find_package(Threads REQUIRED)
target_link_libraries(TCRtrie PRIVATE Threads::Threads CLI11::CLI11)
