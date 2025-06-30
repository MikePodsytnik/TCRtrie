#pragma once

#include <string>
#include <vector>

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
