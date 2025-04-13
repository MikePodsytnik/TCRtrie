#pragma once

#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <cstring>
#include <future>
#include <algorithm>
#include <optional>
#include <fstream>
#include <sstream>
#include "AirrParser.h"

class Trie {
public:
    int insertionCost = 4;
    int depletionCost = 4;

    struct TrieNode {
        std::array<TrieNode*, 30> children{};
        std::vector<int> patternsIndices;
    };

    Trie(const std::string& dataPath);
    Trie();
    Trie(const Trie& other);
    Trie& operator=(const Trie& other);
    Trie(Trie&& other) noexcept;
    Trie& operator=(Trie&& other) noexcept;
    ~Trie();

    std::vector<AIRREntity> Search(const std::string& query, int maxEdits,
                                   const std::optional<std::string>& vGeneFilter = std::nullopt,
                                   const std::optional<std::string>& jGeneFilter = std::nullopt);

    std::vector<AIRREntity> SearchWithScore(const std::string& query, int maxScore,
                                            const std::optional<std::string>& vGeneFilter = std::nullopt,
                                            const std::optional<std::string>& jGeneFilter = std::nullopt);


    bool SearchAny(const std::string& query, int maxEdits);

    std::unordered_map<std::string, std::vector<AIRREntity>> SearchForAIRR(const std::string& pathToAIRR);

    std::unordered_map<std::string, std::vector<AIRREntity>> SearchForAll(const std::vector<std::string>& queries,
                                                                          int maxEdits,
                                                                          const std::optional<std::string>& vGeneFilter = std::nullopt,
                                                                          const std::optional<std::string>& jGeneFilter = std::nullopt);

    std::unordered_map<std::string, std::vector<AIRREntity>> SearchForAllWithScore(const std::vector<std::string>& queries,
                                                                                   int maxScore,
                                                                                   const std::optional<std::string>& vGeneFilter = std::nullopt,
                                                                                   const std::optional<std::string>& jGeneFilter = std::nullopt);

    void LoadSubstitutionMatrix(const std::string& matrixPath);

    void SetMaxQueryLength(int newMaxQueryLength);

private:
    bool useSubstitutionMatrix_ = false;
    int maxQueryLength_ = 32;

    std::unordered_map<char, std::unordered_map<char, float>> substitutionMatrix_;
    TrieNode* root_;

    std::vector<AIRREntity> patterns_;

    void DeleteTrie(TrieNode* node);

    TrieNode* CopyTrie(const TrieNode* node);

    void SearchRecursive(const std::string &query, int maxEdits, const std::string &currentPrefix,
                         TrieNode* node, const int* prevRow, int queryLength,
                         std::vector<AIRREntity>& results,
                         const std::optional<std::string>& vGeneFilter,
                         const std::optional<std::string>& jGeneFilter);

    void SearchRecursiveScore(const std::string &query, int maxScore, const std::string &currentPrefix,
                              TrieNode* node, const int* prevRow, int queryLength,
                              std::vector<AIRREntity>& results,
                              const std::optional<std::string>& vGeneFilter,
                              const std::optional<std::string>& jGeneFilter);

    bool SearchAnyRecursive(const std::string &query, int maxEdits,
                            TrieNode* node, const int* prevRow, int queryLength);

    void LoadAIRRAndBuildTrie(const std::string& dataPath);

    bool IsStandardAminoAcidSequence(const std::string& seq) const;
};
