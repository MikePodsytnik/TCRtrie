#pragma once

#include <string>
#include <utility>
#include <vector>
#include <unordered_map>

struct AIRREntity {
    std::string junctionAA;
    std::string vGene;
    std::string jGene;
    double distance;

    AIRREntity() {}

    AIRREntity(std::string junctionAA_,
               std::string vGene_,
               std::string jGene_,
               double distance_) :
               junctionAA(std::move(junctionAA_)),
               vGene(std::move(vGene_)),
               jGene(std::move(jGene_)),
               distance(distance_){}
};

std::vector<AIRREntity> ParseAIRR(const std::string& filepath);
