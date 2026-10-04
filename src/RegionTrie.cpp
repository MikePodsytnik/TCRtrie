#include "RegionTrie.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
    constexpr float kInf = std::numeric_limits<float>::infinity();
    constexpr float kEps = 1e-6f;

    int OpPriority(AlignmentOpType op) {
        if (op == AlignmentOpType::Match || op == AlignmentOpType::Substitution) return 0;
        if (op == AlignmentOpType::Deletion) return 1;
        return 2;
    }

    struct Cell {
        float cost = kInf;
        int prevI = -1;
        int prevJ = -1;
        AlignmentOpType op = AlignmentOpType::Match;
        bool reachable = false;
    };

    bool Better(float newCost, AlignmentOpType newOp, const Cell& current) {
        if (!current.reachable) return true;
        if (newCost + kEps < current.cost) return true;
        return std::fabs(newCost - current.cost) <= kEps && OpPriority(newOp) < OpPriority(current.op);
    }
}

RegionTrie::RegionTrie(const std::vector<std::string>& sequences,
                       const std::vector<std::vector<Region>>& regions,
                       const std::vector<std::string>& vGenes,
                       const std::vector<std::string>& jGenes)
    : root_(std::make_unique<Node>()), sequences_(sequences), regions_(regions), vGenes_(vGenes), jGenes_(jGenes) {
    if (sequences_.size() != regions_.size()) {
        throw std::runtime_error("RegionTrie received different number of sequences and region markups");
    }
    BuildTrie();
}

void RegionTrie::LoadRegionMatrices(const std::string& vMatrixPath,
                                    const std::string& ndnMatrixPath,
                                    const std::string& jMatrixPath) {
    matrices_[RegionIndex(Region::V)].Load(vMatrixPath);
    matrices_[RegionIndex(Region::NDN)].Load(ndnMatrixPath);
    matrices_[RegionIndex(Region::J)].Load(jMatrixPath);
}

int RegionTrie::ChildIndex(char letter, Region region) {
    if (letter < 'A' || letter > 'Z') throw std::runtime_error("Invalid amino acid letter");
    return RegionIndex(region) * kAlphabetSize + (letter - 'A');
}

const CostMatrix& RegionTrie::MatrixFor(Region region) const {
    return matrices_[RegionIndex(region)];
}

void RegionTrie::BuildTrie() {
    for (int index = 0; index < static_cast<int>(sequences_.size()); ++index) {
        if (sequences_[index].size() != regions_[index].size()) {
            throw std::runtime_error("CDR3 length and region markup length do not match");
        }
        Node* node = root_.get();
        for (std::size_t pos = 0; pos < sequences_[index].size(); ++pos) {
            char letter = sequences_[index][pos];
            Region region = regions_[index][pos];
            int childIndex = ChildIndex(letter, region);
            if (!node->children[childIndex]) {
                node->children[childIndex] = std::make_unique<Node>();
                node->children[childIndex]->letter = letter;
                node->children[childIndex]->region = region;
            }
            node = node->children[childIndex].get();
        }
        node->indices.push_back(index);
    }
}

std::vector<std::pair<std::size_t, float>> RegionTrie::SearchIndices(
    const std::string& query,
    const std::vector<Region>& queryRegions,
    float maxCost,
    const std::optional<std::string>& vGeneFilter,
    const std::optional<std::string>& jGeneFilter) const {

    if (query.size() != queryRegions.size()) {
        throw std::runtime_error("Query CDR3 length and region markup length do not match");
    }

    for (const auto& matrix : matrices_) {
        if (!matrix.Loaded()) throw std::runtime_error("Region matrices are not loaded");
    }

    std::vector<float> initialRow(query.size() + 1, 0.0f);
    for (std::size_t i = 1; i <= query.size(); ++i) {
        Region region = queryRegions[i - 1];
        initialRow[i] = initialRow[i - 1] + MatrixFor(region).Cost(query[i - 1], '-');
    }

    std::vector<std::pair<std::size_t, float>> results;
    SearchRecursive(*root_, query, queryRegions, maxCost, initialRow, results, vGeneFilter, jGeneFilter);
    return results;
}

void RegionTrie::SearchRecursive(const Node& node,
                                 const std::string& query,
                                 const std::vector<Region>& queryRegions,
                                 float maxCost,
                                 const std::vector<float>& prevRow,
                                 std::vector<std::pair<std::size_t, float>>& results,
                                 const std::optional<std::string>& vGeneFilter,
                                 const std::optional<std::string>& jGeneFilter) const {
    if (!node.indices.empty() && prevRow.back() <= maxCost) {
        for (int index : node.indices) {
            if (vGeneFilter && vGenes_[index] != *vGeneFilter) continue;
            if (jGeneFilter && jGenes_[index] != *jGeneFilter) continue;
            results.emplace_back(static_cast<std::size_t>(index), prevRow.back());
        }
    }

    for (const auto& childPtr : node.children) {
        if (!childPtr) continue;
        const Node& child = *childPtr;
        const CostMatrix& targetMatrix = MatrixFor(child.region);
        std::vector<float> nextRow(query.size() + 1, 0.0f);
        nextRow[0] = prevRow[0] + targetMatrix.Cost('-', child.letter);
        float minValue = nextRow[0];

        for (std::size_t j = 1; j <= query.size(); ++j) {
            char queryLetter = query[j - 1];
            Region queryRegion = queryRegions[j - 1];
            const CostMatrix& queryMatrix = MatrixFor(queryRegion);

            float consumeTarget = prevRow[j] + targetMatrix.Cost('-', child.letter);
            float consumeQuery = nextRow[j - 1] + queryMatrix.Cost(queryLetter, '-');
            float consumeBoth = kInf;

            if (queryRegion == child.region) {
                consumeBoth = prevRow[j - 1] + queryMatrix.Cost(queryLetter, child.letter);
            }

            nextRow[j] = std::min({consumeTarget, consumeQuery, consumeBoth});
            minValue = std::min(minValue, nextRow[j]);
        }

        if (minValue <= maxCost) {
            SearchRecursive(child, query, queryRegions, maxCost, nextRow, results, vGeneFilter, jGeneFilter);
        }
    }
}

std::optional<AlignmentResult> RegionTrie::AlignIndexHit(
    const std::string& query,
    const std::vector<Region>& queryRegions,
    std::size_t targetIndex,
    std::optional<float> maxCost) const {

    if (targetIndex >= sequences_.size()) return std::nullopt;
    if (query.size() != queryRegions.size()) throw std::runtime_error("Query CDR3 length and region markup length do not match");

    const std::string& target = sequences_[targetIndex];
    const std::vector<Region>& targetRegions = regions_[targetIndex];
    int qLen = static_cast<int>(query.size());
    int tLen = static_cast<int>(target.size());

    std::vector<std::vector<Cell>> dp(qLen + 1, std::vector<Cell>(tLen + 1));
    dp[0][0].reachable = true;
    dp[0][0].cost = 0.0f;

    for (int i = 0; i <= qLen; ++i) {
        for (int j = 0; j <= tLen; ++j) {
            if (!dp[i][j].reachable) continue;

            if (i < qLen) {
                Region region = queryRegions[i];
                float nextCost = dp[i][j].cost + MatrixFor(region).Cost(query[i], '-');
                if (Better(nextCost, AlignmentOpType::Deletion, dp[i + 1][j])) {
                    dp[i + 1][j] = {nextCost, i, j, AlignmentOpType::Deletion, true};
                }
            }

            if (j < tLen) {
                Region region = targetRegions[j];
                float nextCost = dp[i][j].cost + MatrixFor(region).Cost('-', target[j]);
                if (Better(nextCost, AlignmentOpType::Insertion, dp[i][j + 1])) {
                    dp[i][j + 1] = {nextCost, i, j, AlignmentOpType::Insertion, true};
                }
            }

            if (i < qLen && j < tLen && queryRegions[i] == targetRegions[j]) {
                Region region = queryRegions[i];
                auto op = query[i] == target[j] ? AlignmentOpType::Match : AlignmentOpType::Substitution;
                float nextCost = dp[i][j].cost + MatrixFor(region).Cost(query[i], target[j]);
                if (Better(nextCost, op, dp[i + 1][j + 1])) {
                    dp[i + 1][j + 1] = {nextCost, i, j, op, true};
                }
            }
        }
    }

    if (!dp[qLen][tLen].reachable) return std::nullopt;
    if (maxCost && dp[qLen][tLen].cost > *maxCost + kEps) return std::nullopt;

    AlignmentResult result;
    result.distance = dp[qLen][tLen].cost;

    std::string queryRev;
    std::string targetRev;
    std::string queryRegionRev;
    std::string targetRegionRev;
    int i = qLen;
    int j = tLen;

    while (i > 0 || j > 0) {
        Cell cell = dp[i][j];

        if (cell.op == AlignmentOpType::Deletion) {
            Region region = queryRegions[i - 1];
            ++result.deletions;
            result.regionCosts[RegionIndex(region)] += MatrixFor(region).Cost(query[i - 1], '-');
            queryRev.push_back(query[i - 1]);
            targetRev.push_back('-');
            queryRegionRev.push_back(RegionCode(region));
            targetRegionRev.push_back('-');
        } else if (cell.op == AlignmentOpType::Insertion) {
            Region region = targetRegions[j - 1];
            ++result.insertions;
            result.regionCosts[RegionIndex(region)] += MatrixFor(region).Cost('-', target[j - 1]);
            queryRev.push_back('-');
            targetRev.push_back(target[j - 1]);
            queryRegionRev.push_back('-');
            targetRegionRev.push_back(RegionCode(region));
        } else {
            Region region = queryRegions[i - 1];
            if (cell.op == AlignmentOpType::Substitution) ++result.substitutions;
            result.regionCosts[RegionIndex(region)] += MatrixFor(region).Cost(query[i - 1], target[j - 1]);
            queryRev.push_back(query[i - 1]);
            targetRev.push_back(target[j - 1]);
            queryRegionRev.push_back(RegionCode(region));
            targetRegionRev.push_back(RegionCode(region));
        }

        i = cell.prevI;
        j = cell.prevJ;
    }

    std::reverse(queryRev.begin(), queryRev.end());
    std::reverse(targetRev.begin(), targetRev.end());
    std::reverse(queryRegionRev.begin(), queryRegionRev.end());
    std::reverse(targetRegionRev.begin(), targetRegionRev.end());

    result.queryAligned = std::move(queryRev);
    result.targetAligned = std::move(targetRev);
    result.queryRegionsAligned = std::move(queryRegionRev);
    result.targetRegionsAligned = std::move(targetRegionRev);
    return result;
}
