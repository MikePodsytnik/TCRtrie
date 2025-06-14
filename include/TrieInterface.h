#pragma once

#include <string>

struct SearchConfig {
    std::string inputPath;
    std::string outputPath;
    std::string query;
    std::string inputQueries;
    int nEdits = 0;
    std::string matrixPath;
    float scoreRadius = 0;
    std::string vGene;
    std::string jGene;
};

void RunSearch(const SearchConfig& config);
