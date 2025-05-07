#pragma once

#include <string>
#include <utility>
#include <vector>
#include <unordered_map>

struct AIRREntity {
    std::string junctionAA;
    std::string vGene;
    std::string jGene;

    AIRREntity() {}

    AIRREntity(std::string junctionAA_,
               std::string vGene_,
               std::string jGene_) :
               junctionAA(std::move(junctionAA_)),
               vGene(std::move(vGene_)),
               jGene(std::move(jGene_)) {}
};

std::vector<AIRREntity> ParseAIRR(const std::string& filepath);
