#pragma once

#include <string>

struct SearchConfig {
    std::string inputPath;
    std::string outputPath;
    std::string query;
    std::string inputQueries;
    int maxSubstitution = -1;
    int maxInsertion = -1;
    int maxDeletion = -1;
    std::string matrixPath;
    float costRadius = -1;
    float deletionScore = -6;
    std::string vGene;
    std::string jGene;
};

void RunSearch(const SearchConfig& config);
