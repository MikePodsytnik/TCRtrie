#pragma once

#include "Region.h"

#include <array>
#include <string>

#undef Match

enum class AlignmentOpType {
    Match,
    Substitution,
    Insertion,
    Deletion
};

struct AlignmentResult {
    std::string queryAligned;
    std::string targetAligned;
    std::string queryRegionsAligned;
    std::string targetRegionsAligned;
    int substitutions = 0;
    int insertions = 0;
    int deletions = 0;
    float distance = 0.0f;
    std::array<float, kRegionCount> regionCosts{};
};
