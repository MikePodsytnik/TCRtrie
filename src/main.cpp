#include <iostream>
#include <fstream>
#include "Trie.h"

void load_patterns(const std::string& filename, std::vector<std::string>& patterns) {
    std::ifstream file(filename);
    std::string line;
    while(std::getline(file, line)) {
        patterns.push_back(line);
    }
}

int main() {
    std::vector<std::string> patterns;
    load_patterns("../VDJdb.txt", patterns);
    Trie trie = Trie(patterns);

    auto start_time = std::chrono::high_resolution_clock::now();
    auto result = trie.SearchForAll(patterns, 3);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;
    cout << "All VDJdb: " << elapsed_time.count() << '\n';

    return 0;
}