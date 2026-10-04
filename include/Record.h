#pragma once

#include "Region.h"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

struct Record {
    std::size_t rowIndex = 0;
    std::string junctionAA;
    std::string vGene;
    std::string jGene;
    std::string epitope;
    std::string species;
    std::string chain;
    std::optional<int> vEnd;
    std::optional<int> jStart;
    std::vector<Region> regions;
    std::vector<std::string> rawFields;
};
