#pragma once

#include "AirrParser.h"

#include <array>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

class Trie {
public:
    struct TrieNode {
        std::array<TrieNode*, 30> children{};
        std::vector<int> indices;
    };

    struct Stat {
        int distance;
        int insertion;
        int deletion;
        int substitution;
    };

    explicit Trie(const std::vector<std::string>& sequences);
    explicit Trie(const std::string& dataPath);
    Trie();
    Trie(const Trie& other);
    Trie& operator=(const Trie& other);
    Trie(Trie&& other) noexcept;
    Trie& operator=(Trie&& other) noexcept;
    ~Trie();

    std::vector<std::string> Search(const std::string& query, int maxEdits);

    std::unordered_map<std::string, std::vector<std::string>> Search(const std::vector<std::string>& queries,
                                                                    int maxEdits);

    std::vector<AIRREntity> SearchAIRR(const std::string& query,
                                       int maxSubstitution,
                                       int maxInsertion,
                                       int maxDeletion,
                                       const std::optional<std::string>& vGeneFilter = std::nullopt,
                                       const std::optional<std::string>& jGeneFilter = std::nullopt);

    std::vector<AIRREntity> SearchWithMatrix(const std::string& query, float maxCost,
                                             const std::optional<std::string>& vGeneFilter = std::nullopt,
                                             const std::optional<std::string>& jGeneFilter = std::nullopt);

    bool SearchAny(const std::string& query, int maxEdits);

    std::unordered_map<std::string, std::vector<AIRREntity>> SearchForAll(const std::vector<std::string>& queries,
                                                                          int maxSubstitution,
                                                                          int maxInsertion,
                                                                          int maxDeletion,
                                                                          const std::optional<std::string>& vGeneFilter = std::nullopt,
                                                                          const std::optional<std::string>& jGeneFilter = std::nullopt);

    std::unordered_map<std::string, std::vector<AIRREntity>> SearchForAllWithMatrix(const std::vector<std::string>& queries,
                                                                                    float maxCost,
                                                                                    const std::optional<std::string>& vGeneFilter = std::nullopt,
                                                                                    const std::optional<std::string>& jGeneFilter = std::nullopt);

    void LoadSubstitutionMatrix(const std::string& matrixPath);

    void SetMaxQueryLength(int newMaxQueryLength);

private:
    bool useSubstitutionMatrix_ = false;
    int maxQueryLength_ = 32;

    std::unordered_map<char, std::unordered_map<char, float>> substitutionMatrix_;
    TrieNode* root_;

    std::vector<std::string> sequences_;
    std::vector<std::string> vGenes_;
    std::vector<std::string> jGenes_;

    void DeleteTrie(TrieNode* node);

    TrieNode* CopyTrie(const TrieNode* node);

    void SearchRecursive(const std::string &query, int maxEdits,
                         const std::string &currentPrefix, TrieNode* node,
                         const int* prevRow, int queryLength,
                         std::vector<std::string>& results);

    void SearchRecursiveAIRR(const std::string &query, int maxEdits,
                         TrieNode* node, const int* prevRow, int queryLength,
                         std::vector<AIRREntity>& results,
                         const std::optional<std::string>& vGeneFilter,
                         const std::optional<std::string>& jGeneFilter);

    void SearchRecursiveCost(const std::string &query, float maxCost,
                             TrieNode* node, const float * prevRow, int queryLength,
                             std::vector<AIRREntity>& results,
                             const std::optional<std::string>& vGeneFilter,
                             const std::optional<std::string>& jGeneFilter);

    bool SearchAnyRecursive(const std::string &query, int maxEdits,
                            TrieNode* node, const int* prevRow, int queryLength);

    std::vector<Stat> PruneStats(const std::vector<Stat>& stats);

    std::vector<Stat> DetailedLevenshteinAll(
            const std::string& s,
            const std::string& t,
            int maxEdits);

    void LoadAIRR(const std::string& dataPath);

    void BuildTrie();
};