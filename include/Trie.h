#pragma once

#include "AirrParser.h"

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

class Trie {
public:
    static constexpr int kAlphabetSize = 26;

    struct TrieNode {
        std::array<TrieNode*, kAlphabetSize> children{};
        std::vector<int> indices;
    };

    struct Stat {
        int distance;
        int insertion;
        int deletion;
        int substitution;
    };

    struct EditState {
        int16_t sub;
        int16_t ins;
        int16_t del;
        int total() const { return sub + ins + del; }
    };

    static constexpr int kMaxPareto = 64;

    struct StatCell {
        EditState data[kMaxPareto];
        int size = 0;

        void clear() { size = 0; }
        bool empty() const { return size == 0; }
        void push_back(EditState s) {
            if (size < kMaxPareto) data[size++] = s;
        }
    };

    explicit Trie();
    explicit Trie(const std::string& dataPath);
    explicit Trie(const std::vector<std::string>& sequences,
                  const std::vector<std::string>& vGenes,
                  const std::vector<std::string>& jGenes);

    Trie(const Trie& other);
    Trie& operator=(const Trie& other);

    Trie(Trie&& other) noexcept;
    Trie& operator=(Trie&& other) noexcept;

    ~Trie();

    void SetMaxQueryLength(int newMaxQueryLength);
    void LoadSubstitutionMatrix(const std::string& matrixPath);
    void SetDeletionScore(float deletionScore);
    void PrintMatrix();

    std::vector<std::string> Search(const std::string& query, int maxEdits);

    std::unordered_map<std::string, std::vector<std::string>> Search(
        const std::vector<std::string>& queries,
        int maxEdits);

    std::vector<AIRREntity> SearchAIRR(
        const std::string& query,
        int maxSubstitution = 0,
        int maxInsertion = 0,
        int maxDeletion = 0,
        std::optional<int> maxEdits = std::nullopt,
        const std::optional<std::string>& vGeneFilter = std::nullopt,
        const std::optional<std::string>& jGeneFilter = std::nullopt);

    std::unordered_map<std::string, std::vector<AIRREntity>> SearchForAll(
        const std::vector<std::string>& queries,
        int maxSubstitution = 0,
        int maxInsertion = 0,
        int maxDeletion = 0,
        std::optional<int> maxEdits = std::nullopt,
        std::optional<std::vector<std::string>> vGeneFilters = std::nullopt,
        std::optional<std::vector<std::string>> jGeneFilters = std::nullopt);

    std::vector<AIRREntity> SearchWithMatrix(
        const std::string& query,
        float maxCost,
        const std::optional<std::string>& vGeneFilter = std::nullopt,
        const std::optional<std::string>& jGeneFilter = std::nullopt);

    std::unordered_map<std::string, std::vector<AIRREntity>> SearchForAllWithMatrix(
        const std::vector<std::string>& queries,
        float maxCost,
        std::optional<std::vector<std::string>> vGeneFilters = std::nullopt,
        std::optional<std::vector<std::string>> jGeneFilters = std::nullopt);

    std::vector<std::pair<size_t, int>> SearchIndices(
        const std::string& query,
        int maxSubstitution = 0,
        int maxInsertion = 0,
        int maxDeletion = 0,
        std::optional<int> maxEdits = std::nullopt,
        const std::optional<std::string>& vGeneFilter = std::nullopt,
        const std::optional<std::string>& jGeneFilter = std::nullopt);

    std::unordered_map<std::string, std::vector<std::pair<size_t, int>>> SearchIndicesForAll(
        const std::vector<std::string>& queries,
        int maxSubstitution = 0,
        int maxInsertion = 0,
        int maxDeletion = 0,
        std::optional<int> maxEdits = std::nullopt,
        std::optional<std::vector<std::string>> vGeneFilters = std::nullopt,
        std::optional<std::vector<std::string>> jGeneFilters = std::nullopt);

    std::vector<std::pair<size_t, float>> SearchIndicesWithMatrix(
        const std::string& query,
        float maxCost,
        const std::optional<std::string>& vGeneFilter = std::nullopt,
        const std::optional<std::string>& jGeneFilter = std::nullopt);

    std::unordered_map<std::string, std::vector<std::pair<size_t, float>>> SearchIndicesForAllWithMatrix(
        const std::vector<std::string>& queries,
        float maxCost,
        std::optional<std::vector<std::string>> vGeneFilters = std::nullopt,
        std::optional<std::vector<std::string>> jGeneFilters = std::nullopt);

    std::unordered_set<AIRREntity> ClusterUsage(
        const std::vector<std::string>& cluster,
        int maxSubstitution = 0,
        int maxInsertion = 0,
        int maxDeletion = 0,
        std::optional<int> maxEdits = std::nullopt,
        std::optional<std::vector<std::string>> vGeneFilters = std::nullopt,
        std::optional<std::vector<std::string>> jGeneFilters = std::nullopt);

    std::unordered_set<AIRREntity> ClusterUsageWithMatrix(
        const std::vector<std::string>& cluster,
        float maxCost,
        std::optional<std::vector<std::string>> vGeneFilters = std::nullopt,
        std::optional<std::vector<std::string>> jGeneFilters = std::nullopt);

    bool SearchAny(const std::string& query, int maxEdits);

private:
    bool useSubstitutionMatrix_ = false;
    int maxQueryLength_ = 64;
    float deletionScore_ = -6.0f;

    std::unordered_map<char, std::unordered_map<char, float>> substitutionMatrix_;

    TrieNode* root_ = nullptr;

    std::vector<std::string> sequences_;
    std::vector<std::string> vGenes_;
    std::vector<std::string> jGenes_;

    void LoadAIRR(const std::string& dataPath);
    void BuildTrie();

    void DeleteTrie(TrieNode* node);
    TrieNode* CopyTrie(const TrieNode* node);

    void UpdateSubstitutionMatrix(float deletionScore);

    void SearchSubstitutionOnly(
        const std::string& query,
        int maxSub,
        TrieNode* node,
        int depth,
        int mismatches,
        int queryLength,
        std::vector<AIRREntity>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter);

    void SearchSubstitutionOnlyIDs(
        const std::string& query,
        int maxSub,
        TrieNode* node,
        int depth,
        int mismatches,
        int queryLength,
        std::vector<std::pair<size_t, int>>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter);

    void SearchRecursive(
        const std::string& query,
        int maxEdits,
        TrieNode* node,
        const int* prevRow,
        int queryLength,
        std::vector<std::string>& results);

    void SearchRecursiveAIRR(
        const std::string& query,
        int maxEdits,
        TrieNode* node,
        const int* prevRow,
        int queryLength,
        std::vector<AIRREntity>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter);

    void SearchRecursiveIDs(
        const std::string& query,
        int maxEdits,
        TrieNode* node,
        const int* prevRow,
        int queryLength,
        std::vector<std::pair<size_t, int>>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter);

    template <typename ResultType, typename EmitFunc>
    void SearchRecursiveDetailed(
        const std::string& query,
        int maxEdits,
        int maxSub,
        int maxIns,
        int maxDel,
        TrieNode* node,
        const int* prevRowSimple,
        const StatCell* prevRowDetailed,
        int queryLength,
        std::vector<ResultType>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter,
        EmitFunc emitFunc);

    static void PrunePareto(StatCell& cell);

    void SearchRecursiveCost(
        const std::string& query,
        float maxCost,
        TrieNode* node,
        const float* prevRow,
        int queryLength,
        std::vector<AIRREntity>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter);

    void SearchRecursiveCostIDs(
        const std::string& query,
        float maxCost,
        TrieNode* node,
        const float* prevRow,
        int queryLength,
        std::vector<std::pair<size_t, float>>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter);

    bool SearchAnyRecursive(
        const std::string& query,
        int maxEdits,
        TrieNode* node,
        const int* prevRow,
        int queryLength);
};