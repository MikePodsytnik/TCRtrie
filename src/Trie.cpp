#include "Trie.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <deque>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <thread>

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
        jGenes_(jGenes) {
    BuildTrie();
}

Trie::Trie() : root_(new TrieNode()) {}

Trie::Trie(const Trie& other)
        : root_(nullptr),
          maxQueryLength_(other.maxQueryLength_),
          useSubstitutionMatrix_(other.useSubstitutionMatrix_),
          substitutionMatrix_(other.substitutionMatrix_),
          sequences_(other.sequences_),
          vGenes_(other.vGenes_),
          jGenes_(other.jGenes_) {
    root_ = CopyTrie(other.root_);
}

Trie::Trie(Trie&& other) noexcept
        : root_(other.root_),
          maxQueryLength_(other.maxQueryLength_),
          useSubstitutionMatrix_(other.useSubstitutionMatrix_),
          substitutionMatrix_(std::move(other.substitutionMatrix_)),
          sequences_(std::move(other.sequences_)),
          vGenes_(std::move(other.vGenes_)),
          jGenes_(std::move(other.jGenes_)) {
    other.root_ = nullptr;
}

Trie& Trie::operator=(const Trie& other) {
    if (this != &other) {
        DeleteTrie(root_);
        maxQueryLength_ = other.maxQueryLength_;
        useSubstitutionMatrix_ = other.useSubstitutionMatrix_;
        substitutionMatrix_ = other.substitutionMatrix_;
        sequences_ = other.sequences_;
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
        maxQueryLength_ = other.maxQueryLength_;
        useSubstitutionMatrix_ = other.useSubstitutionMatrix_;
        substitutionMatrix_ = std::move(other.substitutionMatrix_);
        sequences_ = std::move(other.sequences_);
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
    if (queryLength > maxQueryLength_) {
        std::cerr << "Query length exceeds maximum allowed length." << std::endl;
        return results;
    }
    static constexpr int MAX_Q = 33;
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

    if (queryLength > maxQueryLength_) {
        std::cerr << query << " :query length exceeds maximum allowed length("
                  << maxQueryLength_ << ")" << std::endl;
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
        static constexpr int MAX_Q = 33;
        int initialRow[MAX_Q];
        for (int i = 0; i <= queryLength; ++i) initialRow[i] = i;
        SearchRecursiveAIRR(query, *maxEdits, root_, initialRow, queryLength,
                            results, vGeneFilter, jGeneFilter);
        return results;
    }

    static constexpr int MAX_Q = 33;

    int initialRowSimple[MAX_Q];
    for (int i = 0; i <= queryLength; ++i) initialRowSimple[i] = i;

    StatCell initialRowDetailed[MAX_Q];
    initialRowDetailed[0].push_back({0, 0, 0});
    for (int j = 1; j <= queryLength; ++j) {
        if (j <= *maxEdits) {
            initialRowDetailed[j].push_back({0, static_cast<int16_t>(j), 0});
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

    if (queryLength > maxQueryLength_) {
        std::cerr << query << " :query length exceeds maximum allowed length("
                  << maxQueryLength_ << ")" << std::endl;
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
        static constexpr int MAX_Q = 33;
        int initialRow[MAX_Q];
        for (int i = 0; i <= queryLength; ++i) initialRow[i] = i;
        SearchRecursiveIDs(query, *maxEdits, root_, initialRow, queryLength,
                           results, vGeneFilter, jGeneFilter);
        return results;
    }

    static constexpr int MAX_Q = 33;

    int initialRowSimple[MAX_Q];
    for (int i = 0; i <= queryLength; ++i) initialRowSimple[i] = i;

    StatCell initialRowDetailed[MAX_Q];
    initialRowDetailed[0].push_back({0, 0, 0});
    for (int j = 1; j <= queryLength; ++j) {
        if (j <= *maxEdits) {
            initialRowDetailed[j].push_back({0, static_cast<int16_t>(j), 0});
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
    if (queryLength > maxQueryLength_) {
        std::cerr << "Query length exceeds maximum allowed length." << std::endl;
        return false;
    }
    static constexpr int MAX_Q = 33;
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

    if (queryLength > maxQueryLength_) {
        std::cerr << "Query length exceeds maximum allowed length." << std::endl;
        return results;
    }

    static constexpr int MAX_Q = 33;
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
    if (queryLength > maxQueryLength_) {
        std::cerr << "Query length exceeds maximum allowed length." << std::endl;
        return results;
    }

    static constexpr int MAX_Q = 33;
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

std::unordered_map<std::string, std::vector<int>> Trie::SearchGroupIdsForAll(
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

    std::unordered_map<std::string, std::vector<int>> result;
    std::vector<std::future<std::pair<std::string, std::vector<int>>>> futures;

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
             maxEdits, vFilter, jFilter, unique]() -> std::pair<std::string, std::vector<int>> {

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
                    return { query, std::move(gids) };
                }

                std::unordered_set<int> s;
                s.reserve(hits.size());
                for (const auto& [idx, dist] : hits) {
                    s.insert(groupIds_[idx]);
                }

                std::vector<int> gids;
                gids.reserve(s.size());
                for (int gid : s) gids.push_back(gid);

                return { query, std::move(gids) };
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

    static constexpr int MAX_Q = 33;

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

    static constexpr int MAX_Q = 33;

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

    static constexpr int MAX_Q = 33;

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

    static constexpr int MAX_Q = 33;

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

    static constexpr int MAX_Q = 33;

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

    static constexpr int MAX_Q = 33;

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

    static constexpr int MAX_Q = 33;
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

void Trie::LoadSubstitutionMatrix(const std::string& matrixPath) {
    std::ifstream file(matrixPath);
    if (!file) { std::cerr << "Cannot open matrix\n"; return; }

    std::vector<char> letters;
    std::string line;
    std::getline(file, line);
    std::istringstream hs(line);
    char letter;
    while (hs >> letter) letters.push_back(letter);

    std::unordered_map<char, std::unordered_map<char, float>> rawScores;
    bool isCostMatrix = true, isDiagonalMatrix = true;

    rawScores['-']['-'] = std::fabs(deletionScore_);

    for (char r : letters) {
        rawScores[r]['-'] = deletionScore_;
        rawScores['-'][r] = deletionScore_;
        file >> letter;
        for (char c : letters) {
            float v; file >> v;
            rawScores[r][c] = v;
            rawScores[c][r] = v;
            if (v < 0) isCostMatrix = false;
            if (r == c && v != 0) isDiagonalMatrix = false;
        }
    }

    deletionScore_ = rawScores['-']['-'];
    letters.push_back('-');
    substitutionMatrix_.clear();

    if (isCostMatrix && isDiagonalMatrix) {
        substitutionMatrix_ = rawScores;
    } else {
        for (char r : letters) {
            for (char c : letters) {
                substitutionMatrix_[r][c] = (rawScores[r][r] + rawScores[c][c]) * 0.5f - rawScores[r][c];
            }
        }
    }

    for (char r : letters)
        for (char c : letters)
            if (substitutionMatrix_[r][c] < 0)
                std::cerr << "Negative cost: " << r << " vs " << c << " = " << substitutionMatrix_[r][c] << std::endl;

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
}

void Trie::SetMaxQueryLength(int n) {
    maxQueryLength_ = n;
}

void Trie::UpdateSubstitutionMatrix(float deletionScore) {
    std::vector<char> keys;
    keys.reserve(substitutionMatrix_.size());
    for (const auto& kv : substitutionMatrix_) keys.push_back(kv.first);
    for (auto c : keys) {
        if (c == '-') continue;
        substitutionMatrix_[c]['-'] -= deletionScore_ * 0.5f;
        substitutionMatrix_['-'][c] -= deletionScore_ * 0.5f;
        substitutionMatrix_[c]['-'] += std::fabs(deletionScore) * 0.5f;
        substitutionMatrix_['-'][c] += std::fabs(deletionScore) * 0.5f;
    }
}

void Trie::SetDeletionScore(float deletionScore) {
    std::cout << "New deletion score: " << deletionScore << std::endl;
    if (useSubstitutionMatrix_) {
        std::cout << "New Substitution-Score Matrix:" << std::endl;
        UpdateSubstitutionMatrix(deletionScore);
        PrintMatrix();
    }
    deletionScore_ = deletionScore;
}