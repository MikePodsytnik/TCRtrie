#pragma once

#include <string>
#include <vector>
#include <unordered_map>

struct AIRREntity {
    std::string junctionAA;
    std::string vGene;
    std::string jGene;
};

std::vector<AIRREntity> ParseAIRR(const std::string& filepath);
