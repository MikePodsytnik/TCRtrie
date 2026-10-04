#pragma once

#include "Alignment.h"
#include "CostMatrix.h"

#include <array>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

class Trie {
public:
    Trie(const std::vector<std::string>& sequences,
         const std::vector<std::string>& vGenes,
         const std::vector<std::string>& jGenes);

    void LoadSubstitutionMatrix(const std::string& matrixPath);

    std::vector<std::pair<std::size_t, int>> SearchIndices(
        const std::string& query,
        int maxSubstitution,
        int maxInsertion,
        int maxDeletion,
        int maxEdits,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter) const;

    std::vector<std::pair<std::size_t, float>> SearchIndicesWithMatrix(
        const std::string& query,
        float maxCost,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter) const;

    std::optional<AlignmentResult> AlignIndexHit(
        const std::string& query,
        std::size_t targetIndex,
        int maxSubstitution,
        int maxInsertion,
        int maxDeletion,
        int maxEdits) const;

    std::optional<AlignmentResult> AlignIndexHitWithMatrix(
        const std::string& query,
        std::size_t targetIndex,
        std::optional<float> maxCost) const;

private:
    struct Node {
        std::array<std::unique_ptr<Node>, 26> children{};
        std::vector<int> indices;
        char letter = 0;
    };

    std::unique_ptr<Node> root_;
    std::vector<std::string> sequences_;
    std::vector<std::string> vGenes_;
    std::vector<std::string> jGenes_;
    CostMatrix matrix_;

    void BuildTrie();
    void SearchMatrixRecursive(const Node& node,
                               const std::string& query,
                               float maxCost,
                               const std::vector<float>& prevRow,
                               std::vector<std::pair<std::size_t, float>>& results,
                               const std::optional<std::string>& vGeneFilter,
                               const std::optional<std::string>& jGeneFilter) const;

    std::optional<int> BoundedEditDistance(const std::string& query,
                                           const std::string& target,
                                           int maxSubstitution,
                                           int maxInsertion,
                                           int maxDeletion,
                                           int maxEdits) const;
};
