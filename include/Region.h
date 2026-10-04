#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

constexpr int kRegionCount = 3;

enum class Region : std::uint8_t {
    V = 0,
    NDN = 1,
    J = 2
};

inline int RegionIndex(Region region) {
    return static_cast<int>(region);
}

inline char RegionCode(Region region) {
    switch (region) {
        case Region::V:
            return 'V';
        case Region::NDN:
            return 'N';
        case Region::J:
            return 'J';
    }
    return '?';
}

inline Region RegionFromIndex(int index) {
    if (index == 0) return Region::V;
    if (index == 1) return Region::NDN;
    if (index == 2) return Region::J;
    throw std::runtime_error("Invalid region index");
}

inline std::string RegionString(const std::vector<Region>& regions) {
    std::string result;
    result.reserve(regions.size());
    for (Region region : regions) {
        result.push_back(RegionCode(region));
    }
    return result;
}
