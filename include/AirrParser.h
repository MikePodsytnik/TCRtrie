#pragma once

#include <string>
#include <utility>
#include <vector>
#include <unordered_map>

struct AIRREntity {
    std::string junctionAA;
    std::string vGene;
    std::string jGene;
    double distance{};

    AIRREntity() = default;
    AIRREntity(std::string_view ja,
               std::string_view v,
               std::string_view j,
               double d)
            : junctionAA(ja)
            , vGene(v)
            , jGene(j)
            , distance(d)
    {}
};

std::vector<AIRREntity> ParseAIRR(const std::string& filepath);
