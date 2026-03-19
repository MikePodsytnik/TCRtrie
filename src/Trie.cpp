#include "Trie.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cctype>
#include <deque>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace {
    constexpr int MAX_Q = 64;
}

Trie::Trie(const std::string& dataPath) {
    root_ = new TrieNode();
    LoadAIRR(dataPath);
    BuildTrie();
}

Trie::Trie(const std::vector<std::string>& sequences,
           const std::vector<std::string>& vGenes,
           const std::vector<std::string>& jGenes) :
        root_(new TrieNode()),
        sequences_(sequences),
        vGenes_(vGenes),
        jGenes_(jGenes),
        groupIds_(sequences.size(), 0) {
    BuildTrie();
}

Trie::Trie() : root_(new TrieNode()) {}

Trie::Trie(const Trie& other)
        : root_(nullptr),
          useSubstitutionMatrix_(other.useSubstitutionMatrix_),
          substitutionMatrix_(other.substitutionMatrix_),
          sequences_(other.sequences_),
          groupIds_(other.groupIds_),
          vGenes_(other.vGenes_),
          jGenes_(other.jGenes_) {
    root_ = CopyTrie(other.root_);
}

Trie::Trie(Trie&& other) noexcept
        : root_(other.root_),
          useSubstitutionMatrix_(other.useSubstitutionMatrix_),
          substitutionMatrix_(std::move(other.substitutionMatrix_)),
          sequences_(std::move(other.sequences_)),
          groupIds_(std::move(other.groupIds_)),
          vGenes_(std::move(other.vGenes_)),
          jGenes_(std::move(other.jGenes_)) {
    other.root_ = nullptr;
}

Trie& Trie::operator=(const Trie& other) {
    if (this != &other) {
        DeleteTrie(root_);
        useSubstitutionMatrix_ = other.useSubstitutionMatrix_;
        substitutionMatrix_ = other.substitutionMatrix_;
        sequences_ = other.sequences_;
        groupIds_ = other.groupIds_;
        vGenes_ = other.vGenes_;
        jGenes_ = other.jGenes_;
        root_ = CopyTrie(other.root_);
    }
    return *this;
}

Trie& Trie::operator=(Trie&& other) noexcept {
    if (this != &other) {
        DeleteTrie(root_);
        root_ = other.root_;
        useSubstitutionMatrix_ = other.useSubstitutionMatrix_;
        substitutionMatrix_ = std::move(other.substitutionMatrix_);
        sequences_ = std::move(other.sequences_);
        groupIds_ = std::move(other.groupIds_);
        vGenes_ = std::move(other.vGenes_);
        jGenes_ = std::move(other.jGenes_);
        other.root_ = nullptr;
    }
    return *this;
}

Trie::~Trie() {
    DeleteTrie(root_);
}

std::vector<std::string> Trie::Search(const std::string& query, int maxEdits) {
    std::vector<std::string> results;
    int queryLength = static_cast<int>(query.size());
    if (queryLength >= MAX_Q) {
        std::cerr << "Query length exceeds maximum allowed length." << std::endl;
        return results;
    }
    
    int initialRow[MAX_Q];
    for (int i = 0; i <= queryLength; ++i) initialRow[i] = i;
    SearchRecursive(query, maxEdits, root_, initialRow, queryLength, results);
    return results;
}

std::vector<AIRREntity> Trie::SearchAIRR(const std::string& query,
                                         int maxSubstitution,
                                         int maxInsertion,
                                         int maxDeletion,
                                         std::optional<int> maxEdits,
                                         const std::optional<std::string>& vGeneFilter,
                                         const std::optional<std::string>& jGeneFilter) {
    if (!maxEdits.has_value() || *maxEdits < 0) {
        maxEdits = maxSubstitution + maxInsertion + maxDeletion;
    }
    std::vector<AIRREntity> results;
    int queryLength = static_cast<int>(query.size());

    if (queryLength >= MAX_Q) {
        std::cerr << query << " :query length exceeds maximum allowed length("
                  << MAX_Q - 1 << ")" << std::endl;
        return results;
    }

    if (maxInsertion == 0 && maxDeletion == 0) {
        SearchSubstitutionOnly(query, maxSubstitution, root_, 0, 0,
                               queryLength, results, vGeneFilter, jGeneFilter);
        return results;
    }

    if (maxSubstitution >= *maxEdits
        && maxInsertion >= *maxEdits
        && maxDeletion >= *maxEdits) {
        int initialRow[MAX_Q];
        for (int i = 0; i <= queryLength; ++i) initialRow[i] = i;
        SearchRecursiveAIRR(query, *maxEdits, root_, initialRow, queryLength,
                            results, vGeneFilter, jGeneFilter);
        return results;
    }

    int initialRowSimple[MAX_Q];
    for (int i = 0; i <= queryLength; ++i) initialRowSimple[i] = i;

    StatCell initialRowDetailed[MAX_Q];
    initialRowDetailed[0].push_back({0, 0, 0});
    for (int j = 1; j <= queryLength; ++j) {
        if (j <= *maxEdits) {
            initialRowDetailed[j].push_back({0, 0, static_cast<int16_t>(j)});
        }
    }

    auto emitAIRR = [this](std::vector<AIRREntity>& res, int index, int dist) {
        res.emplace_back(sequences_[index], vGenes_[index], jGenes_[index], dist);
    };

    SearchRecursiveDetailed(query, *maxEdits, maxSubstitution, maxInsertion, maxDeletion,
                            root_, initialRowSimple, initialRowDetailed,
                            queryLength, results, vGeneFilter, jGeneFilter, emitAIRR);

    return results;
}

std::vector<std::pair<size_t, int>> Trie::SearchIndices(const std::string& query,
                                                        int maxSubstitution,
                                                        int maxInsertion,
                                                        int maxDeletion,
                                                        std::optional<int> maxEdits,
                                                        const std::optional<std::string>& vGeneFilter,
                                                        const std::optional<std::string>& jGeneFilter) {
    if (!maxEdits.has_value() || *maxEdits < 0) {
        maxEdits = maxSubstitution + maxInsertion + maxDeletion;
    }
    std::vector<std::pair<size_t, int>> results;
    int queryLength = static_cast<int>(query.size());

    if (queryLength >= MAX_Q) {
        std::cerr << query << " :query length exceeds maximum allowed length("
                  << MAX_Q - 1 << ")" << std::endl;
        return results;
    }

    if (maxInsertion == 0 && maxDeletion == 0) {
        SearchSubstitutionOnlyIDs(query, maxSubstitution, root_, 0, 0,
                                  queryLength, results, vGeneFilter, jGeneFilter);
        return results;
    }

    if (maxSubstitution >= *maxEdits
        && maxInsertion >= *maxEdits
        && maxDeletion >= *maxEdits) {
        int initialRow[MAX_Q];
        for (int i = 0; i <= queryLength; ++i) initialRow[i] = i;
        SearchRecursiveIDs(query, *maxEdits, root_, initialRow, queryLength,
                           results, vGeneFilter, jGeneFilter);
        return results;
    }

    int initialRowSimple[MAX_Q];
    for (int i = 0; i <= queryLength; ++i) initialRowSimple[i] = i;

    StatCell initialRowDetailed[MAX_Q];
    initialRowDetailed[0].push_back({0, 0, 0});
    for (int j = 1; j <= queryLength; ++j) {
        if (j <= *maxEdits) {
            initialRowDetailed[j].push_back({0, 0, static_cast<int16_t>(j)});
        }
    }

    auto emitIndex = [](std::vector<std::pair<size_t, int>>& res, int index, int dist) {
        res.emplace_back(index, dist);
    };

    SearchRecursiveDetailed(query, *maxEdits, maxSubstitution, maxInsertion, maxDeletion,
                            root_, initialRowSimple, initialRowDetailed,
                            queryLength, results, vGeneFilter, jGeneFilter, emitIndex);

    return results;
}

bool Trie::SearchAny(const std::string& query, int maxEdits) {
    int queryLength = static_cast<int>(query.size());

    if (queryLength >= MAX_Q) {
        std::cerr << "Query length exceeds maximum allowed length." << std::endl;
        return false;
    }
    int initialRow[MAX_Q];
    for (int i = 0; i <= queryLength; ++i) initialRow[i] = i;
    return SearchAnyRecursive(query, maxEdits, root_, initialRow, queryLength);
}

std::vector<AIRREntity> Trie::SearchWithMatrix(const std::string& query, float maxCost,
                                               const std::optional<std::string>& vGeneFilter,
                                               const std::optional<std::string>& jGeneFilter) {
    std::vector<AIRREntity> results;
    int queryLength = static_cast<int>(query.size());

    if (!useSubstitutionMatrix_) {
        std::cerr << "No substitution matrix is entered, only Levenshtein distance search is available" << std::endl;
        return results;
    }

    if (queryLength >= MAX_Q) {
        std::cerr << "Query length exceeds maximum allowed length." << std::endl;
        return results;
    }

    float initialRow[MAX_Q];
    initialRow[0] = 0.0f;
    for (int i = 1; i <= queryLength; ++i) {
        initialRow[i] = initialRow[i - 1] + substitutionMatrix_.at('-').at(query[i - 1]);
    }

    SearchRecursiveCost(query, maxCost, root_, initialRow, queryLength,
                        results, vGeneFilter, jGeneFilter);

    return results;
}

std::vector<std::pair<size_t, float>> Trie::SearchIndicesWithMatrix(
        const std::string& query, float maxCost,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter) {

    std::vector<std::pair<size_t, float>> results;
    int queryLength = static_cast<int>(query.size());

    if (!useSubstitutionMatrix_) {
        std::cerr << "No substitution matrix is entered, only Levenshtein distance search is available" << std::endl;
        return results;
    }

    if (queryLength >= MAX_Q) {
        std::cerr << "Query length exceeds maximum allowed length." << std::endl;
        return results;
    }

    float initialRow[MAX_Q];
    initialRow[0] = 0.0f;
    for (int i = 1; i <= queryLength; ++i) {
        initialRow[i] = initialRow[i - 1] + substitutionMatrix_.at('-').at(query[i - 1]);
    }

    SearchRecursiveCostIDs(query, maxCost, root_, initialRow, queryLength,
                           results, vGeneFilter, jGeneFilter);

    return results;
}

std::unordered_map<std::string, std::vector<std::string>> Trie::Search(
        const std::vector<std::string>& queries, int maxEdits) {

    std::unordered_map<std::string, std::vector<std::string>> result;
    std::vector<std::future<std::pair<std::string, std::vector<std::string>>>> futures;

    std::size_t hc = std::thread::hardware_concurrency();
    std::size_t maxConcurrent = 10 * (hc == 0 ? 1 : hc);

    for (std::size_t i = 0; i < queries.size(); ++i) {
        const std::string query = queries[i];
        futures.emplace_back(std::async(std::launch::async,
            [this, query, maxEdits]() -> std::pair<std::string, std::vector<std::string>> {
                return { query, this->Search(query, maxEdits) };
            }));

        if (futures.size() >= maxConcurrent || i + 1 == queries.size()) {
            for (auto& fut : futures) {
                auto completed = fut.get();
                result[std::move(completed.first)] = std::move(completed.second);
            }
            futures.clear();
        }
    }

    return result;
}

std::unordered_map<std::string, std::vector<AIRREntity>> Trie::SearchForAll(
        const std::vector<std::string>& queries,
        int maxSubstitution, int maxInsertion,
        int maxDeletion, std::optional<int> maxEdits,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters) {

    if (vGeneFilters && vGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as queries");
    }
    if (jGeneFilters && jGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as queries");
    }

    std::unordered_map<std::string, std::vector<AIRREntity>> result;
    std::vector<std::future<std::pair<std::string, std::vector<AIRREntity>>>> futures;

    std::size_t hc = std::thread::hardware_concurrency();
    std::size_t maxConcurrent = 10 * (hc == 0 ? 1 : hc);

    for (std::size_t i = 0; i < queries.size(); ++i) {
        const std::string query = queries[i];

        const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
        const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;

        futures.emplace_back(std::async(std::launch::async,
            [this, query,
             maxSubstitution, maxInsertion, maxDeletion,
             maxEdits, vFilter, jFilter]() -> std::pair<std::string, std::vector<AIRREntity>> {
                return { query, this->SearchAIRR(query, maxSubstitution, maxInsertion,
                                                maxDeletion, maxEdits, vFilter, jFilter) };
            }));

        if (futures.size() >= maxConcurrent || i + 1 == queries.size()) {
            for (auto& fut : futures) {
                auto completed = fut.get();
                result[std::move(completed.first)] = std::move(completed.second);
            }
            futures.clear();
        }
    }

    return result;
}

std::unordered_map<std::string, std::vector<std::pair<size_t, int>>> Trie::SearchIndicesForAll(
        const std::vector<std::string>& queries,
        int maxSubstitution, int maxInsertion,
        int maxDeletion, std::optional<int> maxEdits,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters) {

    if (vGeneFilters && vGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as queries");
    }
    if (jGeneFilters && jGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as queries");
    }

    std::unordered_map<std::string, std::vector<std::pair<size_t, int>>> result;
    std::vector<std::future<std::pair<std::string, std::vector<std::pair<size_t, int>>>>> futures;

    std::size_t hc = std::thread::hardware_concurrency();
    std::size_t maxConcurrent = 10 * (hc == 0 ? 1 : hc);

    for (std::size_t i = 0; i < queries.size(); ++i) {
        const std::string query = queries[i];

        const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
        const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;

        futures.emplace_back(std::async(std::launch::async,
            [this, query,
             maxSubstitution, maxInsertion, maxDeletion,
             maxEdits, vFilter, jFilter]() -> std::pair<std::string, std::vector<std::pair<size_t, int>>> {
                return { query, this->SearchIndices(query,
                                                    maxSubstitution,
                                                    maxInsertion,
                                                    maxDeletion,
                                                    maxEdits,
                                                    vFilter,
                                                    jFilter) };
            }));

        if (futures.size() >= maxConcurrent || i + 1 == queries.size()) {
            for (auto& fut : futures) {
                auto completed = fut.get();
                result[std::move(completed.first)] = std::move(completed.second);
            }
            futures.clear();
        }
    }

    return result;
}

std::unordered_map<std::string, std::vector<AIRREntity>> Trie::SearchForAllWithMatrix(
        const std::vector<std::string>& queries, float maxCost,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters) {

    if (vGeneFilters && vGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as queries");
    }
    if (jGeneFilters && jGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as queries");
    }

    std::unordered_map<std::string, std::vector<AIRREntity>> result;
    std::vector<std::future<std::pair<std::string, std::vector<AIRREntity>>>> futures;

    std::size_t hc = std::thread::hardware_concurrency();
    std::size_t maxConcurrent = 10 * (hc == 0 ? 1 : hc);

    for (std::size_t i = 0; i < queries.size(); ++i) {
        const std::string query = queries[i];

        const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
        const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;

        futures.emplace_back(std::async(std::launch::async,
            [this, query, maxCost, vFilter, jFilter]() -> std::pair<std::string, std::vector<AIRREntity>> {
                return { query, this->SearchWithMatrix(query, maxCost, vFilter, jFilter) };
            }));

        if (futures.size() >= maxConcurrent || i + 1 == queries.size()) {
            for (auto& fut : futures) {
                auto completed = fut.get();
                result[std::move(completed.first)] = std::move(completed.second);
            }
            futures.clear();
        }
    }

    return result;
}

std::unordered_map<std::string, std::vector<std::pair<size_t, float>>> Trie::SearchIndicesForAllWithMatrix(
        const std::vector<std::string>& queries, float maxCost,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters) {

    if (vGeneFilters && vGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as queries");
    }
    if (jGeneFilters && jGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as queries");
    }

    std::unordered_map<std::string, std::vector<std::pair<size_t, float>>> result;
    std::vector<std::future<std::pair<std::string, std::vector<std::pair<size_t, float>>>>> futures;

    std::size_t hc = std::thread::hardware_concurrency();
    std::size_t maxConcurrent = 10 * (hc == 0 ? 1 : hc);

    for (std::size_t i = 0; i < queries.size(); ++i) {
        const std::string query = queries[i];

        const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
        const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;

        futures.emplace_back(std::async(std::launch::async,
            [this, query, maxCost, vFilter, jFilter]()
                -> std::pair<std::string, std::vector<std::pair<size_t, float>>> {
                return { query, this->SearchIndicesWithMatrix(query, maxCost, vFilter, jFilter) };
            }));

        if (futures.size() >= maxConcurrent || i + 1 == queries.size()) {
            for (auto& fut : futures) {
                auto completed = fut.get();
                result[std::move(completed.first)] = std::move(completed.second);
            }
            futures.clear();
        }
    }

    return result;
}

std::vector<std::vector<int>> Trie::SearchGroupIdsForAll(
    const std::vector<std::string>& queries,
    int maxSubstitution, int maxInsertion,
    int maxDeletion, std::optional<int> maxEdits,
    std::optional<std::vector<std::string>> vGeneFilters,
    std::optional<std::vector<std::string>> jGeneFilters,
    bool unique) {
    if (vGeneFilters && vGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as queries");
    }
    if (jGeneFilters && jGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as queries");
    }

    const std::size_t n = queries.size();
    std::vector<std::vector<int>> result(n);

    std::vector<std::future<void>> futures;

    std::size_t hc = std::thread::hardware_concurrency();
    std::size_t maxConcurrent = 10 * (hc == 0 ? 1 : hc);

    for (std::size_t i = 0; i < n; ++i) {

        const std::string query = queries[i];

        const std::optional<std::string> vFilter =
            vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;

        const std::optional<std::string> jFilter =
            jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;

        futures.emplace_back(std::async(std::launch::async,
            [this, &result, i, query,
             maxSubstitution, maxInsertion, maxDeletion,
             maxEdits, vFilter, jFilter, unique]() {

                auto hits = this->SearchIndices(query,
                                                maxSubstitution,
                                                maxInsertion,
                                                maxDeletion,
                                                maxEdits,
                                                vFilter,
                                                jFilter);

                if (!unique) {
                    std::vector<int> gids;
                    gids.reserve(hits.size());
                    for (const auto& [idx, dist] : hits) {
                        gids.push_back(groupIds_[idx]);
                    }
                    result[i] = std::move(gids);
                    return;
                }

                std::unordered_set<int> s;
                s.reserve(hits.size());
                for (const auto& [idx, dist] : hits) {
                    s.insert(groupIds_[idx]);
                }

                std::vector<int> gids;
                gids.reserve(s.size());
                for (int gid : s) {
                    gids.push_back(gid);
                }

                result[i] = std::move(gids);
            }));

        if (futures.size() >= maxConcurrent || i + 1 == n) {
            for (auto& fut : futures) {
                fut.get();
            }
            futures.clear();
        }
    }

    return result;
}

std::unordered_set<AIRREntity> Trie::ClusterUsage(
        const std::vector<std::string>& cluster,
        int maxSubstitution, int maxInsertion, int maxDeletion,
        std::optional<int> maxEdits,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters) {

    if (vGeneFilters && vGeneFilters->size() != cluster.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as cluster");
    }
    if (jGeneFilters && jGeneFilters->size() != cluster.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as cluster");
    }

    std::unordered_set<AIRREntity> result;
    std::vector<std::future<std::vector<AIRREntity>>> futures;

    std::size_t hc = std::thread::hardware_concurrency();
    std::size_t maxConcurrent = 10 * (hc == 0 ? 1 : hc);

    for (std::size_t i = 0; i < cluster.size(); ++i) {
        const std::string query = cluster[i];

        const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
        const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;

        futures.emplace_back(std::async(std::launch::async,
            [this, query,
             maxSubstitution, maxInsertion, maxDeletion,
             maxEdits, vFilter, jFilter]() -> std::vector<AIRREntity> {
                return this->SearchAIRR(query, maxSubstitution, maxInsertion,
                                        maxDeletion, maxEdits, vFilter, jFilter);
            }));

        if (futures.size() >= maxConcurrent) {
            for (auto& fut : futures) {
                auto matches = fut.get();
                for (auto& entity : matches) {
                    result.insert(std::move(entity));
                }
            }
            futures.clear();
        }
    }

    for (auto& fut : futures) {
        auto matches = fut.get();
        for (auto& entity : matches) {
            result.insert(std::move(entity));
        }
    }

    return result;
}

std::unordered_set<AIRREntity> Trie::ClusterUsageWithMatrix(
        const std::vector<std::string>& cluster,
        float maxCost,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters) {

    if (vGeneFilters && vGeneFilters->size() != cluster.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as cluster");
    }
    if (jGeneFilters && jGeneFilters->size() != cluster.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as cluster");
    }

    std::unordered_set<AIRREntity> result;
    std::vector<std::future<std::vector<AIRREntity>>> futures;

    std::size_t hc = std::thread::hardware_concurrency();
    std::size_t maxConcurrent = 10 * (hc == 0 ? 1 : hc);

    for (std::size_t i = 0; i < cluster.size(); ++i) {
        const std::string query = cluster[i];

        const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
        const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;

        futures.emplace_back(std::async(std::launch::async,
            [this, query, maxCost, vFilter, jFilter]() -> std::vector<AIRREntity> {
                return this->SearchWithMatrix(query, maxCost, vFilter, jFilter);
            }));

        if (futures.size() >= maxConcurrent) {
            for (auto& fut : futures) {
                auto matches = fut.get();
                for (auto& entity : matches) {
                    result.insert(std::move(entity));
                }
            }
            futures.clear();
        }
    }

    for (auto& fut : futures) {
        auto matches = fut.get();
        for (auto& entity : matches) {
            result.insert(std::move(entity));
        }
    }

    return result;
}

void Trie::PrunePareto(StatCell& cell) {
    StatCell result;
    for (int i = 0; i < cell.size; ++i) {
        const auto& s = cell.data[i];
        bool dominated = false;
        for (int k = 0; k < result.size; ++k) {
            const auto& r = result.data[k];
            if (r.sub <= s.sub && r.ins <= s.ins && r.del <= s.del) {
                dominated = true;
                break;
            }
        }
        if (dominated) continue;
        StatCell filtered;
        for (int k = 0; k < result.size; ++k) {
            const auto& r = result.data[k];
            if (!(s.sub <= r.sub && s.ins <= r.ins && s.del <= r.del)) {
                filtered.push_back(r);
            }
        }
        filtered.push_back(s);
        result = filtered;
    }
    cell = result;
}

void Trie::SearchSubstitutionOnly(const std::string& query, int maxSub,
                                  TrieNode* node, int depth, int mismatches,
                                  int queryLength,
                                  std::vector<AIRREntity>& results,
                                  const std::optional<std::string>& vGeneFilter,
                                  const std::optional<std::string>& jGeneFilter) {
    if (depth == queryLength) {
        if (!node->indices.empty()) {
            for (int index : node->indices) {
                bool vMatch = !vGeneFilter || vGenes_[index] == *vGeneFilter;
                bool jMatch = !jGeneFilter || jGenes_[index] == *jGeneFilter;
                if (vMatch && jMatch) {
                    results.emplace_back(sequences_[index],
                                         vGenes_[index],
                                         jGenes_[index],
                                         mismatches);
                }
            }
        }
        return;
    }

    char queryChar = query[depth];
    for (int ci = 0; ci < kAlphabetSize; ++ci) {
        TrieNode* child = node->children[ci];
        if (!child) continue;

        int newMismatches = mismatches + (('A' + ci) != queryChar ? 1 : 0);
        if (newMismatches > maxSub) continue;

        SearchSubstitutionOnly(query, maxSub, child, depth + 1, newMismatches,
                               queryLength, results, vGeneFilter, jGeneFilter);
    }
}

void Trie::SearchSubstitutionOnlyIDs(const std::string& query, int maxSub,
                                     TrieNode* node, int depth, int mismatches,
                                     int queryLength,
                                     std::vector<std::pair<size_t, int>>& results,
                                     const std::optional<std::string>& vGeneFilter,
                                     const std::optional<std::string>& jGeneFilter) {
    if (depth == queryLength) {
        if (!node->indices.empty()) {
            for (int index : node->indices) {
                bool vMatch = !vGeneFilter || vGenes_[index] == *vGeneFilter;
                bool jMatch = !jGeneFilter || jGenes_[index] == *jGeneFilter;
                if (vMatch && jMatch) {
                    results.emplace_back(static_cast<size_t>(index), mismatches);
                }
            }
        }
        return;
    }

    char queryChar = query[depth];
    for (int ci = 0; ci < kAlphabetSize; ++ci) {
        TrieNode* child = node->children[ci];
        if (!child) continue;

        int newMismatches = mismatches + (('A' + ci) != queryChar ? 1 : 0);
        if (newMismatches > maxSub) continue;

        SearchSubstitutionOnlyIDs(query, maxSub, child, depth + 1, newMismatches,
                                  queryLength, results, vGeneFilter, jGeneFilter);
    }
}

void Trie::SearchRecursive(const std::string& query, int maxEdits,
                           TrieNode* node, const int* prevRow, int queryLength,
                           std::vector<std::string>& results) {
    if (!node->indices.empty() && prevRow[queryLength] <= maxEdits) {
        for (int index : node->indices) {
            results.push_back(sequences_[index]);
        }
    }

    int minVal = prevRow[0];
    for (int j = 1; j <= queryLength; ++j) {
        if (prevRow[j] < minVal) minVal = prevRow[j];
    }
    if (minVal > maxEdits) return;

    for (int ci = 0; ci < kAlphabetSize; ++ci) {
        TrieNode* child = node->children[ci];
        if (!child) continue;
        char letter = static_cast<char>('A' + ci);

        int nextRow[MAX_Q];
        nextRow[0] = prevRow[0] + 1;
        for (int j = 1; j <= queryLength; ++j) {
            int cost = (query[j - 1] == letter) ? 0 : 1;
            int d = prevRow[j] + 1;
            int i = nextRow[j - 1] + 1;
            int s = prevRow[j - 1] + cost;
            int val = d;
            if (i < val) val = i;
            if (s < val) val = s;
            nextRow[j] = val;
        }
        SearchRecursive(query, maxEdits, child, nextRow, queryLength, results);
    }
}

void Trie::SearchRecursiveAIRR(const std::string& query, int maxEdits,
                               TrieNode* node, const int* prevRow, int queryLength,
                               std::vector<AIRREntity>& results,
                               const std::optional<std::string>& vGeneFilter,
                               const std::optional<std::string>& jGeneFilter) {
    if (!node->indices.empty() && prevRow[queryLength] <= maxEdits) {
        for (int index : node->indices) {
            bool vMatch = !vGeneFilter || vGenes_[index] == *vGeneFilter;
            bool jMatch = !jGeneFilter || jGenes_[index] == *jGeneFilter;
            if (vMatch && jMatch) {
                results.emplace_back(sequences_[index],
                                     vGenes_[index],
                                     jGenes_[index],
                                     prevRow[queryLength]);
            }
        }
    }

    int minVal = prevRow[0];
    for (int j = 1; j <= queryLength; ++j) {
        if (prevRow[j] < minVal) minVal = prevRow[j];
    }
    if (minVal > maxEdits) return;

    for (int ci = 0; ci < kAlphabetSize; ++ci) {
        TrieNode* child = node->children[ci];
        if (!child) continue;
        char letter = static_cast<char>('A' + ci);

        int nextRow[MAX_Q];
        nextRow[0] = prevRow[0] + 1;
        for (int j = 1; j <= queryLength; ++j) {
            int cost = (query[j - 1] == letter) ? 0 : 1;
            int d = prevRow[j] + 1;
            int i = nextRow[j - 1] + 1;
            int s = prevRow[j - 1] + cost;
            int val = d;
            if (i < val) val = i;
            if (s < val) val = s;
            nextRow[j] = val;
        }

        SearchRecursiveAIRR(query, maxEdits, child, nextRow, queryLength,
                            results, vGeneFilter, jGeneFilter);
    }
}

void Trie::SearchRecursiveIDs(const std::string& query, int maxEdits,
                              TrieNode* node, const int* prevRow, int queryLength,
                              std::vector<std::pair<size_t, int>>& results,
                              const std::optional<std::string>& vGeneFilter,
                              const std::optional<std::string>& jGeneFilter) {
    if (!node->indices.empty() && prevRow[queryLength] <= maxEdits) {
        for (int index : node->indices) {
            bool vMatch = !vGeneFilter || vGenes_[index] == *vGeneFilter;
            bool jMatch = !jGeneFilter || jGenes_[index] == *jGeneFilter;
            if (vMatch && jMatch) {
                results.emplace_back(static_cast<size_t>(index), prevRow[queryLength]);
            }
        }
    }

    int minVal = prevRow[0];
    for (int j = 1; j <= queryLength; ++j) {
        if (prevRow[j] < minVal) minVal = prevRow[j];
    }
    if (minVal > maxEdits) return;

    for (int ci = 0; ci < kAlphabetSize; ++ci) {
        TrieNode* child = node->children[ci];
        if (!child) continue;
        char letter = static_cast<char>('A' + ci);

        int nextRow[MAX_Q];
        nextRow[0] = prevRow[0] + 1;
        for (int j = 1; j <= queryLength; ++j) {
            int cost = (query[j - 1] == letter) ? 0 : 1;
            int d = prevRow[j] + 1;
            int i = nextRow[j - 1] + 1;
            int s = prevRow[j - 1] + cost;
            int val = d;
            if (i < val) val = i;
            if (s < val) val = s;
            nextRow[j] = val;
        }

        SearchRecursiveIDs(query, maxEdits, child, nextRow, queryLength,
                           results, vGeneFilter, jGeneFilter);
    }
}

template<typename ResultType, typename EmitFunc>
void Trie::SearchRecursiveDetailed(
        const std::string& query,
        int maxEdits, int maxSub, int maxIns, int maxDel,
        TrieNode* node,
        const int* prevRowSimple,
        const StatCell* prevRowDetailed,
        int queryLength,
        std::vector<ResultType>& results,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter,
        EmitFunc emitFunc) {
    if (!node->indices.empty() && prevRowSimple[queryLength] <= maxEdits) {
        bool found = false;
        const auto& cell = prevRowDetailed[queryLength];
        for (int i = 0; i < cell.size; ++i) {
            const auto& st = cell.data[i];
            if (st.sub <= maxSub && st.ins <= maxIns && st.del <= maxDel) {
                found = true;
                break;
            }
        }
        if (found) {
            for (int index : node->indices) {
                bool vMatch = !vGeneFilter || vGenes_[index] == *vGeneFilter;
                bool jMatch = !jGeneFilter || jGenes_[index] == *jGeneFilter;
                if (vMatch && jMatch) {
                    emitFunc(results, index, prevRowSimple[queryLength]);
                }
            }
        }
    }

    int minVal = prevRowSimple[0];
    for (int j = 1; j <= queryLength; ++j) {
        if (prevRowSimple[j] < minVal) minVal = prevRowSimple[j];
    }
    if (minVal > maxEdits) return;

    for (int ci = 0; ci < kAlphabetSize; ++ci) {
        TrieNode* child = node->children[ci];
        if (!child) continue;
        char letter = static_cast<char>('A' + ci);

        int nextRowSimple[MAX_Q];
        nextRowSimple[0] = prevRowSimple[0] + 1;
        int nextMinVal = nextRowSimple[0];

        for (int j = 1; j <= queryLength; ++j) {
            int cost = (query[j - 1] == letter) ? 0 : 1;
            int d = prevRowSimple[j] + 1;
            int i = nextRowSimple[j - 1] + 1;
            int s = prevRowSimple[j - 1] + cost;
            int val = d;
            if (i < val) val = i;
            if (s < val) val = s;
            nextRowSimple[j] = val;
            if (val < nextMinVal) nextMinVal = val;
        }

        if (nextMinVal > maxEdits) continue;

        StatCell nextRowDetailed[MAX_Q];

        {
            const auto& pcell = prevRowDetailed[0];
            for (int k = 0; k < pcell.size; ++k) {
                const auto& st = pcell.data[k];
                EditState ns = {st.sub, static_cast<int16_t>(st.ins + 1), st.del};
                if (ns.total() <= maxEdits) {
                    nextRowDetailed[0].push_back(ns);
                }
            }
            PrunePareto(nextRowDetailed[0]);
        }

        for (int j = 1; j <= queryLength; ++j) {
            if (nextRowSimple[j] > maxEdits) continue;

            StatCell& cand = nextRowDetailed[j];

            {
                const auto& pcell = prevRowDetailed[j];
                for (int k = 0; k < pcell.size; ++k) {
                    const auto& st = pcell.data[k];
                    EditState ns = {st.sub, static_cast<int16_t>(st.ins + 1), st.del};
                    if (ns.total() <= maxEdits)
                        cand.push_back(ns);
                }
            }

            {
                const auto& ncell = nextRowDetailed[j - 1];
                for (int k = 0; k < ncell.size; ++k) {
                    const auto& st = ncell.data[k];
                    EditState ns = {st.sub, st.ins, static_cast<int16_t>(st.del + 1)};
                    if (ns.total() <= maxEdits)
                        cand.push_back(ns);
                }
            }

            {
                int cost = (query[j - 1] == letter) ? 0 : 1;
                const auto& pcell = prevRowDetailed[j - 1];
                for (int k = 0; k < pcell.size; ++k) {
                    const auto& st = pcell.data[k];
                    EditState ns = {static_cast<int16_t>(st.sub + cost), st.ins, st.del};
                    if (ns.total() <= maxEdits)
                        cand.push_back(ns);
                }
            }

            PrunePareto(cand);
        }

        SearchRecursiveDetailed(query, maxEdits, maxSub, maxIns, maxDel,
                                child, nextRowSimple, nextRowDetailed,
                                queryLength, results, vGeneFilter, jGeneFilter,
                                emitFunc);
    }
}

void Trie::SearchRecursiveCost(const std::string& query, float maxCost,
                               TrieNode* node, const float* prevRow, int queryLength,
                               std::vector<AIRREntity>& results,
                               const std::optional<std::string>& vGeneFilter,
                               const std::optional<std::string>& jGeneFilter) {
    if (!node->indices.empty() && (prevRow[queryLength] <= maxCost)) {
        for (int index : node->indices) {
            bool vMatch = !vGeneFilter || vGenes_[index] == *vGeneFilter;
            bool jMatch = !jGeneFilter || jGenes_[index] == *jGeneFilter;
            if (vMatch && jMatch) {
                results.emplace_back(sequences_[index],
                                     vGenes_[index],
                                     jGenes_[index],
                                     prevRow[queryLength]);
            }
        }
    }

    for (int ci = 0; ci < kAlphabetSize; ++ci) {
        TrieNode* child = node->children[ci];
        if (!child) continue;
        char letter = static_cast<char>('A' + ci);

        float nextRow[MAX_Q];
        float deletionCost = substitutionMatrix_.at('-').at(letter);
        nextRow[0] = prevRow[0] + deletionCost;
        float minVal = nextRow[0];

        for (int j = 1; j <= queryLength; ++j) {
            char queryChar = query[j - 1];
            float subCost = substitutionMatrix_.at(queryChar).at(letter);
            float insertionCost = substitutionMatrix_.at('-').at(queryChar);

            float d = prevRow[j] + deletionCost;
            float i = nextRow[j - 1] + insertionCost;
            float s = prevRow[j - 1] + subCost;
            float val = d;
            if (i < val) val = i;
            if (s < val) val = s;
            nextRow[j] = val;
            if (val < minVal) minVal = val;
        }

        if (minVal > maxCost) continue;

        SearchRecursiveCost(query, maxCost, child, nextRow, queryLength,
                            results, vGeneFilter, jGeneFilter);
    }
}

void Trie::SearchRecursiveCostIDs(const std::string& query, float maxCost,
                                  TrieNode* node, const float* prevRow, int queryLength,
                                  std::vector<std::pair<size_t, float>>& results,
                                  const std::optional<std::string>& vGeneFilter,
                                  const std::optional<std::string>& jGeneFilter) {
    if (!node->indices.empty() && (prevRow[queryLength] <= maxCost)) {
        for (int index : node->indices) {
            bool vMatch = !vGeneFilter || vGenes_[index] == *vGeneFilter;
            bool jMatch = !jGeneFilter || jGenes_[index] == *jGeneFilter;
            if (vMatch && jMatch) {
                results.emplace_back(static_cast<size_t>(index), prevRow[queryLength]);
            }
        }
    }

    for (int ci = 0; ci < kAlphabetSize; ++ci) {
        TrieNode* child = node->children[ci];
        if (!child) continue;

        char letter = static_cast<char>('A' + ci);

        float nextRow[MAX_Q];
        float deletionCost = substitutionMatrix_.at('-').at(letter);

        nextRow[0] = prevRow[0] + deletionCost;
        float minVal = nextRow[0];

        for (int j = 1; j <= queryLength; ++j) {
            char queryChar = query[j - 1];

            float subCost = substitutionMatrix_.at(queryChar).at(letter);
            float insertionCost = substitutionMatrix_.at('-').at(queryChar);

            float d = prevRow[j] + deletionCost;
            float i = nextRow[j - 1] + insertionCost;
            float s = prevRow[j - 1] + subCost;

            float val = d;
            if (i < val) val = i;
            if (s < val) val = s;

            nextRow[j] = val;
            if (val < minVal) minVal = val;
        }

        if (minVal > maxCost) continue;

        SearchRecursiveCostIDs(query, maxCost, child, nextRow, queryLength,
                               results, vGeneFilter, jGeneFilter);
    }
}

bool Trie::SearchAnyRecursive(const std::string& query, int maxEdits,
                              TrieNode* node, const int* prevRow, int queryLength) {
    if (!node->indices.empty() && prevRow[queryLength] <= maxEdits) {
        return true;
    }
    int minVal = prevRow[0];
    for (int j = 1; j <= queryLength; ++j) {
        if (prevRow[j] < minVal) minVal = prevRow[j];
    }
    if (minVal > maxEdits) return false;

    for (int ci = 0; ci < kAlphabetSize; ++ci) {
        TrieNode* child = node->children[ci];
        if (!child) continue;
        char letter = static_cast<char>('A' + ci);
        int nextRow[MAX_Q];
        nextRow[0] = prevRow[0] + 1;
        for (int j = 1; j <= queryLength; ++j) {
            int cost = (query[j - 1] == letter) ? 0 : 1;
            int d = prevRow[j] + 1;
            int i = nextRow[j - 1] + 1;
            int s = prevRow[j - 1] + cost;
            int val = d;
            if (i < val) val = i;
            if (s < val) val = s;
            nextRow[j] = val;
        }
        if (SearchAnyRecursive(query, maxEdits, child, nextRow, queryLength))
            return true;
    }
    return false;
}

void Trie::LoadAIRR(const std::string& dataPath) {
    auto entries = ParseAIRR(dataPath);
    for (auto& e : entries) {
        sequences_.push_back(e.junctionAA);
        vGenes_.push_back(e.vGene);
        jGenes_.push_back(e.jGene);
        groupIds_.push_back(e.groupId);
    }
}

void Trie::BuildTrie() {
    for (int idx = 0; idx < static_cast<int>(sequences_.size()); ++idx) {
        const auto& seq = sequences_[idx];
        TrieNode* node = root_;
        for (char c : seq) {
            if (c < 'A' || c > 'Z') continue;
            int i = c - 'A';
            if (!node->children[i]) {
                node->children[i] = new TrieNode();
            }
            node = node->children[i];
        }
        node->indices.push_back(idx);
    }
}

void Trie::DeleteTrie(TrieNode* node) {
    if (!node) return;
    for (auto* child : node->children) {
        DeleteTrie(child);
    }
    delete node;
}

Trie::TrieNode* Trie::CopyTrie(const TrieNode* node) {
    if (!node) return nullptr;
    TrieNode* newNode = new TrieNode();
    newNode->indices = node->indices;
    for (int i = 0; i < kAlphabetSize; ++i) {
        if (node->children[i]) {
            newNode->children[i] = CopyTrie(node->children[i]);
        }
    }
    return newNode;
}

namespace {
    const std::string kAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    std::string Trim(const std::string& s) {
        std::size_t l = 0;
        while (l < s.size() && std::isspace(static_cast<unsigned char>(s[l]))) {
            ++l;
        }

        std::size_t r = s.size();
        while (r > l && std::isspace(static_cast<unsigned char>(s[r - 1]))) {
            --r;
        }

        return s.substr(l, r - l);
    }

    void RemoveBom(std::string& s) {
        if (s.size() >= 3 &&
            static_cast<unsigned char>(s[0]) == 0xEF &&
            static_cast<unsigned char>(s[1]) == 0xBB &&
            static_cast<unsigned char>(s[2]) == 0xBF) {
            s.erase(0, 3);
        }
    }

    std::vector<std::string> SplitLine(const std::string& line, const std::string& delimiter) {
        std::vector<std::string> result;

        if (delimiter.empty()) {
            std::istringstream iss(line);
            std::string token;
            while (iss >> token) {
                result.push_back(token);
            }
            return result;
        }

        std::size_t start = 0;
        while (true) {
            std::size_t pos = line.find(delimiter, start);
            if (pos == std::string::npos) {
                result.push_back(Trim(line.substr(start)));
                break;
            }
            result.push_back(Trim(line.substr(start, pos - start)));
            start = pos + delimiter.size();
        }

        return result;
    }

    char ParseLabel(std::string token) {
        if (token.size() >= 3 &&
            static_cast<unsigned char>(token[0]) == 0xEF &&
            static_cast<unsigned char>(token[1]) == 0xBB &&
            static_cast<unsigned char>(token[2]) == 0xBF) {
            token.erase(0, 3);
        }

        token = Trim(token);

        if (token.size() != 1) {
            throw std::runtime_error("Invalid matrix label: '" + token + "'");
        }

        return static_cast<char>(std::toupper(static_cast<unsigned char>(token[0])));
    }

    std::string Join(const std::vector<char>& values) {
        std::ostringstream out;
        for (std::size_t i = 0; i < values.size(); ++i) {
            if (i) {
                out << ", ";
            }
            out << values[i];
        }
        return out.str();
    }
}

void Trie::LoadSubstitutionMatrix(const std::string& matrixPath,
                                  const std::string& delimiter,
                                  float gapFactor) {
    if (gapFactor < 1.0f) {
        std::cerr << "gapFactor must be >= 1.0\n";
        throw std::runtime_error("Invalid gapFactor");
    }

    std::ifstream file(matrixPath);
    if (!file) {
        std::cerr << "Cannot open matrix: " << matrixPath << "\n";
        throw std::runtime_error("Cannot open matrix");
    }

    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Matrix file is empty\n";
        throw std::runtime_error("Empty matrix file");
    }

    std::vector<std::string> header = SplitLine(line, delimiter);
    if (header.empty()) {
        std::cerr << "Matrix header is empty\n";
        throw std::runtime_error("Invalid matrix header");
    }
    RemoveBom(header[0]);

    if (delimiter.empty()) {
        if (header.size() != 20 && header.size() != 21) {
            std::cerr << "Matrix must contain 20 or 21 columns, got " << header.size() << "\n";
            throw std::runtime_error("Invalid matrix size");
        }
    } else {
        if (header.empty() || !Trim(header[0]).empty()) {
            std::cerr << "Top-left matrix cell must be empty\n";
            throw std::runtime_error("Invalid matrix header");
        }
        header.erase(header.begin());

        if (header.size() != 20 && header.size() != 21) {
            std::cerr << "Matrix must contain 20 or 21 columns, got " << header.size() << "\n";
            throw std::runtime_error("Invalid matrix size");
        }
    }

    std::vector<char> cols;
    for (const std::string& token : header) {
        cols.push_back(ParseLabel(token));
    }

    std::set<char> colSet(cols.begin(), cols.end());
    if (colSet.size() != cols.size()) {
        std::cerr << "Duplicate column labels found\n";
        throw std::runtime_error("Duplicate column labels");
    }

    for (char c : cols) {
        if (kAminoAcids.find(c) == std::string::npos && c != '-') {
            std::cerr << "Unexpected column label: " << c << "\n";
            throw std::runtime_error("Unexpected column label");
        }
    }

    std::vector<char> missing;
    for (char aa : kAminoAcids) {
        if (colSet.find(aa) == colSet.end()) {
            missing.push_back(aa);
        }
    }

    if (!missing.empty()) {
        std::cerr << "Missing amino acids: " << Join(missing) << "\n";
        throw std::runtime_error("Missing amino acids in matrix");
    }

    std::unordered_map<char, std::unordered_map<char, float>> score;
    std::vector<char> rows;

    while (std::getline(file, line)) {
        if (Trim(line).empty()) {
            continue;
        }

        std::vector<std::string> parts = SplitLine(line, delimiter);
        if (parts.size() != cols.size() + 1) {
            std::cerr << "Invalid row width\n";
            throw std::runtime_error("Invalid matrix row");
        }

        char row = ParseLabel(parts[0]);
        rows.push_back(row);

        for (std::size_t i = 0; i < cols.size(); ++i) {
            try {
                std::string token = Trim(parts[i + 1]);
                std::size_t pos = 0;
                float value = std::stof(token, &pos);
                if (pos != token.size() || !std::isfinite(value)) {
                    throw std::runtime_error("");
                }
                score[row][cols[i]] = value;
            } catch (...) {
                std::cerr << "Invalid numeric value in row " << row << ", column " << cols[i] << "\n";
                throw std::runtime_error("Invalid numeric value in matrix");
            }
        }
    }

    if (rows.size() != cols.size()) {
        std::cerr << "Matrix is not square\n";
        throw std::runtime_error("Matrix is not square");
    }

    std::set<char> rowSet(rows.begin(), rows.end());
    if (rowSet.size() != rows.size()) {
        std::cerr << "Duplicate row labels found\n";
        throw std::runtime_error("Duplicate row labels");
    }

    for (char r : rows) {
        if (kAminoAcids.find(r) == std::string::npos && r != '-') {
            std::cerr << "Unexpected row label: " << r << "\n";
            throw std::runtime_error("Unexpected row label");
        }
    }

    missing.clear();
    for (char aa : kAminoAcids) {
        if (rowSet.find(aa) == rowSet.end()) {
            missing.push_back(aa);
        }
    }

    if (!missing.empty()) {
        std::cerr << "Missing amino acids: " << Join(missing) << "\n";
        throw std::runtime_error("Missing amino acids in matrix");
    }

    if (rowSet != colSet) {
        std::cerr << "Row and column labels do not match\n";
        throw std::runtime_error("Row and column labels do not match");
    }

    bool hasGap = colSet.find('-') != colSet.end();

    auto validateDiagonalDominance = [&](const std::vector<char>& labels) {
        std::vector<char> bad;

        for (char aa : kAminoAcids) {
            float diag = score.at(aa).at(aa);
            bool ok = true;

            for (char other : labels) {
                if (other == aa) {
                    continue;
                }

                if (diag <= score.at(aa).at(other) || diag <= score.at(other).at(aa)) {
                    ok = false;
                    break;
                }
            }

            if (!ok) {
                bad.push_back(aa);
            }
        }

        if (!bad.empty()) {
            std::cerr << "Diagonal is not strictly greater for: " << Join(bad) << "\n";
            throw std::runtime_error("Invalid diagonal values");
        }
    };

    validateDiagonalDominance(cols);

    std::vector<char> labels = cols;

    if (!hasGap) {
        score['-']['-'] = 0.0f;

        for (char aa : kAminoAcids) {
            float minValue = std::numeric_limits<float>::infinity();

            for (char rowAa : kAminoAcids) {
                minValue = std::min(minValue, score.at(rowAa).at(aa));
            }

            float gapScore = (minValue < 0.0f)
                                 ? minValue * gapFactor
                                 : minValue / gapFactor;

            score[aa]['-'] = gapScore;
            score['-'][aa] = gapScore;
        }

        labels.push_back('-');
        validateDiagonalDominance(labels);
    } else {
        if (score.find('-') == score.end() ||
            score.at('-').find('-') == score.at('-').end()) {
            std::cerr << "Gap row is incomplete\n";
            throw std::runtime_error("Incomplete gap row");
        }

        for (char aa : labels) {
            if (score.at('-').find(aa) == score.at('-').end() ||
                score.at(aa).find('-') == score.at(aa).end()) {
                std::cerr << "Gap row or column is incomplete for amino acid " << aa << "\n";
                throw std::runtime_error("Incomplete gap row/column");
            }
        }
    }

    substitutionMatrix_.clear();

    substitutionMatrix_.clear();

    for (char r : labels) {
        for (char c : labels) {
            float cost = 0.0f;

            if (r == '-' && c == '-') {
                cost = 0.0f;
            } else if (r == '-') {
                cost = score.at(c).at(c) + std::abs(score.at(r).at(c));
            } else if (c == '-') {
                cost = score.at(r).at(r) + std::abs(score.at(r).at(c));
            } else {
                cost = (score.at(r).at(r) + score.at(c).at(c)) * 0.5f - score.at(r).at(c);
            }

            if (cost < 0.0f) {
                std::cerr << "Negative cost after conversion for pair " << r << ", " << c << "\n";
                throw std::runtime_error("Negative cost after conversion");
            }

            substitutionMatrix_[r][c] = cost;
        }
    }

    useSubstitutionMatrix_ = true;
}

void Trie::PrintMatrix() {
    std::vector<char> keys;
    keys.reserve(substitutionMatrix_.size());
    for (const auto& kv : substitutionMatrix_) keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());
    std::cout << std::setw(4) << "";
    for (char c : keys) std::cout << std::setw(6) << c;
    std::cout << "\n";
    for (char row : keys) {
        std::cout << std::setw(4) << row;
        for (char col : keys)
            std::cout << std::setw(6) << std::fixed << std::setprecision(2) << substitutionMatrix_[row][col];
        std::cout << "\n";
    }
    std::cout << std::endl;
}
