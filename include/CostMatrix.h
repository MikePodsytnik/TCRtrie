#pragma once

#include <string>
#include <unordered_map>

class CostMatrix {
public:
    void Load(const std::string& path, const std::string& delimiter = "", float gapFactor = 1.5f);
    float Cost(char row, char column) const;
    bool Loaded() const;

private:
    std::unordered_map<char, std::unordered_map<char, float>> cost_;
};
