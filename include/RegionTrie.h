#pragma once

#include "Alignment.h"
#include "CostMatrix.h"
#include "Region.h"

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

class RegionTrie {
public:
    RegionTrie(const std::vector<std::string>& sequences,
               const std::vector<std::vector<Region>>& regions,
               const std::vector<std::string>& vGenes,
               const std::vector<std::string>& jGenes);

    void LoadRegionMatrices(const std::string& vMatrixPath,
                            const std::string& ndnMatrixPath,
                            const std::string& jMatrixPath);

    std::vector<std::pair<std::size_t, float>> SearchIndices(
        const std::string& query,
        const std::vector<Region>& queryRegions,
        float maxCost,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter) const;

    std::optional<AlignmentResult> AlignIndexHit(
        const std::string& query,
        const std::vector<Region>& queryRegions,
        std::size_t targetIndex,
        std::optional<float> maxCost) const;

private:
    static constexpr int kAlphabetSize = 26;
    static constexpr int kChildCount = kAlphabetSize * kRegionCount;

    struct Node {
        std::array<std::unique_ptr<Node>, kChildCount> children{};
        std::vector<int> indices;
        char letter = 0;
        Region region = Region::NDN;
    };

    std::unique_ptr<Node> root_;
    std::vector<std::string> sequences_;
    std::vector<std::vector<Region>> regions_;
    std::vector<std::string> vGenes_;
    std::vector<std::string> jGenes_;
    std::array<CostMatrix, kRegionCount> matrices_;

    void BuildTrie();
    static int ChildIndex(char letter, Region region);
    const CostMatrix& MatrixFor(Region region) const;

    void SearchRecursive(const Node& node,
                         const std::string& query,
                         const std::vector<Region>& queryRegions,
                         float maxCost,
                         const std::vector<float>& prevRow,
                         std::vector<std::pair<std::size_t, float>>& results,
                         const std::optional<std::string>& vGeneFilter,
                         const std::optional<std::string>& jGeneFilter) const;
};
