#include "Trie.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <cctype>
#include <deque>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>

namespace {
    constexpr int MAX_Q = 64;
    constexpr float kFloatEps = 1e-6f;

    const std::string kAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    bool IsCanonicalAminoAcidSequence(const std::string& sequence) {
        if (sequence.empty()) {
            return false;
        }

        for (char c : sequence) {
            if (kAminoAcids.find(c) == std::string::npos) {
                return false;
            }
        }

        return true;
    }

    std::size_t ResolveThreadCount(std::optional<std::size_t> numThreads, std::size_t taskCount = 0) {
        std::size_t resolved = numThreads.value_or(4);
        if (resolved == 0) {
            resolved = 1;
        }
        if (taskCount != 0 && resolved > taskCount) {
            resolved = taskCount;
        }
        return resolved;
    }

    template <typename Output, typename Worker>
    std::vector<Output> RunParallelOrdered(std::size_t taskCount,
                                           std::optional<std::size_t> numThreads,
                                           Worker&& worker) {
        if (taskCount == 0) {
            return {};
        }

        const std::size_t threadsCount = ResolveThreadCount(numThreads, taskCount);
        std::vector<Output> result(taskCount);
        std::atomic<std::size_t> next{0};
        std::exception_ptr firstException;
        std::mutex exceptionMutex;
        std::vector<std::thread> threads;
        threads.reserve(threadsCount);

        for (std::size_t t = 0; t < threadsCount; ++t) {
            threads.emplace_back([&]() {
                try {
                    while (true) {
                        const std::size_t i = next.fetch_add(1);
                        if (i >= taskCount) {
                            break;
                        }
                        result[i] = worker(i);
                    }
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lock(exceptionMutex);
                        if (!firstException) {
                            firstException = std::current_exception();
                        }
                    }
                    next.store(taskCount);
                }
            });
        }

        for (auto& thread : threads) {
            thread.join();
        }

        if (firstException) {
            std::rethrow_exception(firstException);
        }

        return result;
    }

    struct BoundedTraceState {
        int16_t sub = 0;
        int16_t ins = 0;
        int16_t del = 0;
        int prevI = -1;
        int prevJ = -1;
        int prevK = -1;
        Trie::AlignmentOpType op = Trie::AlignmentOpType::Match;

        int total() const {
            return static_cast<int>(sub) + static_cast<int>(ins) + static_cast<int>(del);
        }
    };

    bool SameCounts(const BoundedTraceState& lhs, const BoundedTraceState& rhs) {
        return lhs.sub == rhs.sub && lhs.ins == rhs.ins && lhs.del == rhs.del;
    }

    bool Dominates(const BoundedTraceState& lhs, const BoundedTraceState& rhs) {
        return lhs.sub <= rhs.sub && lhs.ins <= rhs.ins && lhs.del <= rhs.del &&
               (lhs.sub < rhs.sub || lhs.ins < rhs.ins || lhs.del < rhs.del);
    }

    void AddPrunedState(std::vector<BoundedTraceState>& cell, const BoundedTraceState& candidate) {
        for (const auto& existing : cell) {
            if (SameCounts(existing, candidate) || Dominates(existing, candidate)) {
                return;
            }
        }

        cell.erase(std::remove_if(cell.begin(), cell.end(), [&](const BoundedTraceState& existing) {
            return Dominates(candidate, existing);
        }), cell.end());

        cell.push_back(candidate);
    }

    struct MatrixTraceCell {
        float cost = std::numeric_limits<float>::infinity();
        int prevI = -1;
        int prevJ = -1;
        Trie::AlignmentOpType op = Trie::AlignmentOpType::Match;
        bool reachable = false;
    };

    int AlignmentOpPriority(Trie::AlignmentOpType op) {
        switch (op) {
            case Trie::AlignmentOpType::Match:
            case Trie::AlignmentOpType::Substitution:
                return 0;
            case Trie::AlignmentOpType::Deletion:
                return 1;
            case Trie::AlignmentOpType::Insertion:
                return 2;
        }
        return 3;
    }

    bool BetterMatrixCandidate(float newCost,
                               Trie::AlignmentOpType newOp,
                               const MatrixTraceCell& current) {
        if (!current.reachable) {
            return true;
        }
        if (newCost + kFloatEps < current.cost) {
            return true;
        }
        if (std::fabs(newCost - current.cost) <= kFloatEps &&
            AlignmentOpPriority(newOp) < AlignmentOpPriority(current.op)) {
            return true;
        }
        return false;
    }
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
    if (sequences.size() != vGenes.size() || sequences.size() != jGenes.size()) {
        throw std::invalid_argument("sequences, vGenes and jGenes must have the same length");
    }

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
        int effectiveMaxSub = std::min(maxSubstitution, *maxEdits);
        SearchSubstitutionOnly(query, effectiveMaxSub, root_, 0, 0,
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
        if (j <= maxDeletion && j <= *maxEdits) {
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
        int effectiveMaxSub = std::min(maxSubstitution, *maxEdits);
        SearchSubstitutionOnlyIDs(query, effectiveMaxSub, root_, 0, 0,
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
        if (j <= maxDeletion && j <= *maxEdits) {
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

    for (char c : query) {
        if (substitutionMatrix_.find(c) == substitutionMatrix_.end()) {
            std::cerr << "Query contains character absent from substitution matrix: "
                      << c << std::endl;
            return results;
        }
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

    for (char c : query) {
        if (substitutionMatrix_.find(c) == substitutionMatrix_.end()) {
            std::cerr << "Query contains character absent from substitution matrix: "
                      << c << std::endl;
            return results;
        }
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

std::vector<std::vector<std::string>> Trie::Search(
        const std::vector<std::string>& queries,
        int maxEdits,
        std::optional<std::size_t> numThreads) {
    return RunParallelOrdered<std::vector<std::string>>(
        queries.size(),
        numThreads,
        [&](std::size_t i) {
            return this->Search(queries[i], maxEdits);
        });
}

std::vector<std::vector<AIRREntity>> Trie::SearchForAll(
        const std::vector<std::string>& queries,
        int maxSubstitution,
        int maxInsertion,
        int maxDeletion,
        std::optional<int> maxEdits,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters,
        std::optional<std::size_t> numThreads) {

    if (vGeneFilters && vGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as queries");
    }
    if (jGeneFilters && jGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as queries");
    }

    return RunParallelOrdered<std::vector<AIRREntity>>(
        queries.size(),
        numThreads,
        [&](std::size_t i) {
            const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
            const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;

            return this->SearchAIRR(
                queries[i],
                maxSubstitution,
                maxInsertion,
                maxDeletion,
                maxEdits,
                vFilter,
                jFilter);
        });
}

std::vector<std::vector<std::pair<size_t, int>>> Trie::SearchIndicesForAll(
        const std::vector<std::string>& queries,
        int maxSubstitution,
        int maxInsertion,
        int maxDeletion,
        std::optional<int> maxEdits,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters,
        std::optional<std::size_t> numThreads) {

    if (vGeneFilters && vGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as queries");
    }
    if (jGeneFilters && jGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as queries");
    }

    return RunParallelOrdered<std::vector<std::pair<size_t, int>>>(
        queries.size(),
        numThreads,
        [&](std::size_t i) {
            const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
            const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;

            return this->SearchIndices(
                queries[i],
                maxSubstitution,
                maxInsertion,
                maxDeletion,
                maxEdits,
                vFilter,
                jFilter);
        });
}

std::vector<std::vector<AIRREntity>> Trie::SearchForAllWithMatrix(
        const std::vector<std::string>& queries,
        float maxCost,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters,
        std::optional<std::size_t> numThreads) {

    if (vGeneFilters && vGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as queries");
    }
    if (jGeneFilters && jGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as queries");
    }

    return RunParallelOrdered<std::vector<AIRREntity>>(
        queries.size(),
        numThreads,
        [&](std::size_t i) {
            const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
            const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;
            return this->SearchWithMatrix(queries[i], maxCost, vFilter, jFilter);
        });
}

std::vector<std::vector<std::pair<size_t, float>>> Trie::SearchIndicesForAllWithMatrix(
        const std::vector<std::string>& queries,
        float maxCost,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters,
        std::optional<std::size_t> numThreads) {

    if (vGeneFilters && vGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as queries");
    }
    if (jGeneFilters && jGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as queries");
    }

    return RunParallelOrdered<std::vector<std::pair<size_t, float>>>(
        queries.size(),
        numThreads,
        [&](std::size_t i) {
            const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
            const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;
            return this->SearchIndicesWithMatrix(queries[i], maxCost, vFilter, jFilter);
        });
}

std::vector<std::vector<int>> Trie::SearchGroupIdsForAll(
    const std::vector<std::string>& queries,
    int maxSubstitution,
    int maxInsertion,
    int maxDeletion,
    std::optional<int> maxEdits,
    std::optional<std::vector<std::string>> vGeneFilters,
    std::optional<std::vector<std::string>> jGeneFilters,
    bool unique,
    std::optional<std::size_t> numThreads) {
    if (vGeneFilters && vGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as queries");
    }
    if (jGeneFilters && jGeneFilters->size() != queries.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as queries");
    }

    return RunParallelOrdered<std::vector<int>>(
        queries.size(),
        numThreads,
        [&](std::size_t i) {
            const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
            const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;

            auto hits = this->SearchIndices(
                queries[i],
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
                return gids;
            }

            std::unordered_set<int> seen;
            seen.reserve(hits.size());
            for (const auto& [idx, dist] : hits) {
                seen.insert(groupIds_[idx]);
            }

            std::vector<int> gids;
            gids.reserve(seen.size());
            for (int gid : seen) {
                gids.push_back(gid);
            }
            return gids;
        });
}

std::unordered_set<AIRREntity> Trie::ClusterUsage(
        const std::vector<std::string>& cluster,
        int maxSubstitution,
        int maxInsertion,
        int maxDeletion,
        std::optional<int> maxEdits,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters,
        std::optional<std::size_t> numThreads) {

    if (vGeneFilters && vGeneFilters->size() != cluster.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as cluster");
    }
    if (jGeneFilters && jGeneFilters->size() != cluster.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as cluster");
    }

    auto batches = RunParallelOrdered<std::vector<AIRREntity>>(
        cluster.size(),
        numThreads,
        [&](std::size_t i) {
            const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
            const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;

            return this->SearchAIRR(
                cluster[i],
                maxSubstitution,
                maxInsertion,
                maxDeletion,
                maxEdits,
                vFilter,
                jFilter);
        });

    std::unordered_set<AIRREntity> result;
    for (auto& batch : batches) {
        for (auto& entity : batch) {
            result.insert(std::move(entity));
        }
    }
    return result;
}

std::unordered_set<AIRREntity> Trie::ClusterUsageWithMatrix(
        const std::vector<std::string>& cluster,
        float maxCost,
        std::optional<std::vector<std::string>> vGeneFilters,
        std::optional<std::vector<std::string>> jGeneFilters,
        std::optional<std::size_t> numThreads) {

    if (vGeneFilters && vGeneFilters->size() != cluster.size()) {
        throw std::invalid_argument("vGeneFilters must have the same length as cluster");
    }
    if (jGeneFilters && jGeneFilters->size() != cluster.size()) {
        throw std::invalid_argument("jGeneFilters must have the same length as cluster");
    }

    auto batches = RunParallelOrdered<std::vector<AIRREntity>>(
        cluster.size(),
        numThreads,
        [&](std::size_t i) {
            const std::optional<std::string> vFilter =
                vGeneFilters ? std::optional<std::string>((*vGeneFilters)[i]) : std::nullopt;
            const std::optional<std::string> jFilter =
                jGeneFilters ? std::optional<std::string>((*jGeneFilters)[i]) : std::nullopt;
            return this->SearchWithMatrix(cluster[i], maxCost, vFilter, jFilter);
        });

    std::unordered_set<AIRREntity> result;
    for (auto& batch : batches) {
        for (auto& entity : batch) {
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
        int bestValidDist = maxEdits + 1;
        const auto& cell = prevRowDetailed[queryLength];

        for (int i = 0; i < cell.size; ++i) {
            const auto& st = cell.data[i];
            if (st.sub <= maxSub && st.ins <= maxIns && st.del <= maxDel) {
                if (st.total() < bestValidDist) {
                    bestValidDist = st.total();
                }
            }
        }

        if (bestValidDist <= maxEdits) {
            for (int index : node->indices) {
                bool vMatch = !vGeneFilter || vGenes_[index] == *vGeneFilter;
                bool jMatch = !jGeneFilter || jGenes_[index] == *jGeneFilter;
                if (vMatch && jMatch) {
                    emitFunc(results, index, bestValidDist);
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

        auto nextRowSimplePtr = std::make_unique<int[]>(MAX_Q);
        int* nextRowSimple = nextRowSimplePtr.get();
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

        auto nextRowDetailedPtr = std::make_unique<StatCell[]>(MAX_Q);
        StatCell* nextRowDetailed = nextRowDetailedPtr.get();

        {
            StatCell& cand = nextRowDetailed[0];
            const auto& pcell = prevRowDetailed[0];
            for (int k = 0; k < pcell.size; ++k) {
                const auto& st = pcell.data[k];
                EditState ns = {st.sub, static_cast<int16_t>(st.ins + 1), st.del};
                if (ns.sub <= maxSub && ns.ins <= maxIns && ns.del <= maxDel && ns.total() <= maxEdits) {
                    cand.push_back(ns);
                }
            }
            PrunePareto(cand);
        }

        for (int j = 1; j <= queryLength; ++j) {
            if (nextRowSimple[j] > maxEdits) continue;

            StatCell& cand = nextRowDetailed[j];

            auto add_cand = [&](EditState s) {
                bool dominated = false;
                for (int k = 0; k < cand.size; ++k) {
                    const auto& r = cand.data[k];
                    if (r.sub <= s.sub && r.ins <= s.ins && r.del <= s.del) {
                        dominated = true;
                        break;
                    }
                }
                if (dominated) return;

                StatCell filtered;
                for (int k = 0; k < cand.size; ++k) {
                    const auto& r = cand.data[k];
                    if (!(s.sub <= r.sub && s.ins <= r.ins && s.del <= r.del)) {
                        filtered.push_back(r);
                    }
                }
                filtered.push_back(s);
                cand = filtered;
            };

            {
                const auto& pcell = prevRowDetailed[j];
                for (int k = 0; k < pcell.size; ++k) {
                    const auto& st = pcell.data[k];
                    EditState ns = {st.sub, static_cast<int16_t>(st.ins + 1), st.del};
                    if (ns.sub <= maxSub && ns.ins <= maxIns && ns.del <= maxDel && ns.total() <= maxEdits)
                        add_cand(ns);
                }
            }

            {
                const auto& ncell = nextRowDetailed[j - 1];
                for (int k = 0; k < ncell.size; ++k) {
                    const auto& st = ncell.data[k];
                    EditState ns = {st.sub, st.ins, static_cast<int16_t>(st.del + 1)};
                    if (ns.sub <= maxSub && ns.ins <= maxIns && ns.del <= maxDel && ns.total() <= maxEdits)
                        add_cand(ns);
                }
            }

            {
                int cost = (query[j - 1] == letter) ? 0 : 1;
                const auto& pcell = prevRowDetailed[j - 1];
                for (int k = 0; k < pcell.size; ++k) {
                    const auto& st = pcell.data[k];
                    EditState ns = {static_cast<int16_t>(st.sub + cost), st.ins, st.del};
                    if (ns.sub <= maxSub && ns.ins <= maxIns && ns.del <= maxDel && ns.total() <= maxEdits)
                        add_cand(ns);
                }
            }
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

        if (seq.size() > MAX_Q) {
            std::cerr << "Skipping too long sequence "
                      << seq << ": length = " << seq.size() << std::endl;
            continue;
        }

        if (!IsCanonicalAminoAcidSequence(seq)) {
            std::cerr << "Skipping sequence with non-canonical amino acid: " << seq << std::endl;
            continue;
        }

        TrieNode* node = root_;
        for (char c : seq) {
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
    if (gapFactor < 0.0f) {
        std::cerr << "gapFactor must be >= 0.0\n";
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
            float gapScore = -score[aa][aa] * gapFactor;

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


std::optional<Trie::AlignmentResult> Trie::AlignQueryToTarget(
        const std::string& query,
        const std::string& target,
        std::optional<int> maxSubstitution,
        std::optional<int> maxInsertion,
        std::optional<int> maxDeletion,
        std::optional<int> maxEdits) {
    const int qLen = static_cast<int>(query.size());
    const int tLen = static_cast<int>(target.size());

    const int maxSubLimit = maxSubstitution.value_or(std::min(qLen, tLen));
    const int maxInsLimit = maxInsertion.value_or(tLen);
    const int maxDelLimit = maxDeletion.value_or(qLen);
    const int maxTotal = maxEdits.value_or(maxSubLimit + maxInsLimit + maxDelLimit);

    if (maxSubLimit < 0 || maxInsLimit < 0 || maxDelLimit < 0 || maxTotal < 0) {
        throw std::invalid_argument("Alignment limits must be non-negative");
    }

    std::vector<std::vector<std::vector<BoundedTraceState>>> cells(
        qLen + 1,
        std::vector<std::vector<BoundedTraceState>>(tLen + 1));

    cells[0][0].push_back({0, 0, 0, -1, -1, -1, AlignmentOpType::Match});

    for (int i = 0; i <= qLen; ++i) {
        for (int j = 0; j <= tLen; ++j) {
            const auto currentStates = cells[i][j];
            for (int k = 0; k < static_cast<int>(currentStates.size()); ++k) {
                const auto& state = currentStates[k];

                if (i < qLen) {
                    BoundedTraceState candidate = state;
                    candidate.del = static_cast<int16_t>(candidate.del + 1);
                    candidate.prevI = i;
                    candidate.prevJ = j;
                    candidate.prevK = k;
                    candidate.op = AlignmentOpType::Deletion;

                    if (candidate.del <= maxDelLimit && candidate.total() <= maxTotal) {
                        AddPrunedState(cells[i + 1][j], candidate);
                    }
                }

                if (j < tLen) {
                    BoundedTraceState candidate = state;
                    candidate.ins = static_cast<int16_t>(candidate.ins + 1);
                    candidate.prevI = i;
                    candidate.prevJ = j;
                    candidate.prevK = k;
                    candidate.op = AlignmentOpType::Insertion;

                    if (candidate.ins <= maxInsLimit && candidate.total() <= maxTotal) {
                        AddPrunedState(cells[i][j + 1], candidate);
                    }
                }

                if (i < qLen && j < tLen) {
                    BoundedTraceState candidate = state;
                    const bool isMatch = query[i] == target[j];
                    if (!isMatch) {
                        candidate.sub = static_cast<int16_t>(candidate.sub + 1);
                    }
                    candidate.prevI = i;
                    candidate.prevJ = j;
                    candidate.prevK = k;
                    candidate.op = isMatch ? AlignmentOpType::Match : AlignmentOpType::Substitution;

                    if (candidate.sub <= maxSubLimit && candidate.total() <= maxTotal) {
                        AddPrunedState(cells[i + 1][j + 1], candidate);
                    }
                }
            }
        }
    }

    const auto& finalStates = cells[qLen][tLen];
    int bestIndex = -1;
    for (int k = 0; k < static_cast<int>(finalStates.size()); ++k) {
        const auto& st = finalStates[k];
        if (st.sub > maxSubLimit || st.ins > maxInsLimit || st.del > maxDelLimit || st.total() > maxTotal) {
            continue;
        }
        if (bestIndex < 0) {
            bestIndex = k;
            continue;
        }
        const auto& best = finalStates[bestIndex];
        const auto candidateKey = std::make_tuple(st.total(), st.ins + st.del, st.sub, st.del, st.ins);
        const auto bestKey = std::make_tuple(best.total(), best.ins + best.del, best.sub, best.del, best.ins);
        if (candidateKey < bestKey) {
            bestIndex = k;
        }
    }

    if (bestIndex < 0) {
        return std::nullopt;
    }

    AlignmentResult result;
    result.substitutions = finalStates[bestIndex].sub;
    result.insertions = finalStates[bestIndex].ins;
    result.deletions = finalStates[bestIndex].del;
    result.distance = static_cast<float>(finalStates[bestIndex].total());

    std::string queryAlignedRev;
    std::string targetAlignedRev;
    std::vector<AlignmentOp> opsRev;

    int i = qLen;
    int j = tLen;
    int k = bestIndex;

    while (true) {
        const auto& state = cells[i][j][k];
        if (state.prevI < 0) {
            break;
        }

        switch (state.op) {
            case AlignmentOpType::Match:
                queryAlignedRev.push_back(query[i - 1]);
                targetAlignedRev.push_back(target[j - 1]);
                break;
            case AlignmentOpType::Substitution:
                queryAlignedRev.push_back(query[i - 1]);
                targetAlignedRev.push_back(target[j - 1]);
                opsRev.push_back({AlignmentOpType::Substitution, i, query[i - 1], target[j - 1]});
                break;
            case AlignmentOpType::Deletion:
                queryAlignedRev.push_back(query[i - 1]);
                targetAlignedRev.push_back('-');
                opsRev.push_back({AlignmentOpType::Deletion, i, query[i - 1], '-'});
                break;
            case AlignmentOpType::Insertion:
                queryAlignedRev.push_back('-');
                targetAlignedRev.push_back(target[j - 1]);
                opsRev.push_back({AlignmentOpType::Insertion, i, '-', target[j - 1]});
                break;
        }

        const int nextI = state.prevI;
        const int nextJ = state.prevJ;
        const int nextK = state.prevK;
        i = nextI;
        j = nextJ;
        k = nextK;
    }

    std::reverse(queryAlignedRev.begin(), queryAlignedRev.end());
    std::reverse(targetAlignedRev.begin(), targetAlignedRev.end());
    std::reverse(opsRev.begin(), opsRev.end());

    result.queryAligned = std::move(queryAlignedRev);
    result.targetAligned = std::move(targetAlignedRev);
    result.ops = std::move(opsRev);
    return result;
}

std::optional<Trie::AlignmentResult> Trie::AlignIndexHit(
        const std::string& query,
        size_t targetIndex,
        std::optional<int> maxSubstitution,
        std::optional<int> maxInsertion,
        std::optional<int> maxDeletion,
        std::optional<int> maxEdits) {
    if (targetIndex >= sequences_.size()) {
        return std::nullopt;
    }
    return AlignQueryToTarget(
        query,
        sequences_[targetIndex],
        maxSubstitution,
        maxInsertion,
        maxDeletion,
        maxEdits);
}

std::vector<std::optional<Trie::AlignmentResult>> Trie::AlignIndexHits(
        const std::string& query,
        const std::vector<std::pair<size_t, int>>& hits,
        std::optional<int> maxSubstitution,
        std::optional<int> maxInsertion,
        std::optional<int> maxDeletion,
        std::optional<int> maxEdits,
        std::optional<std::size_t> numThreads) {
    return RunParallelOrdered<std::optional<AlignmentResult>>(
        hits.size(),
        numThreads,
        [&](std::size_t i) {
            const auto& [targetIndex, dist] = hits[i];
            (void)dist;
            return AlignIndexHit(
                query,
                targetIndex,
                maxSubstitution,
                maxInsertion,
                maxDeletion,
                maxEdits);
        });
}

std::optional<Trie::AlignmentResult> Trie::AlignQueryToTargetWithMatrix(
        const std::string& query,
        const std::string& target,
        std::optional<float> maxCost) {
    if (!useSubstitutionMatrix_) {
        throw std::runtime_error("No substitution matrix is loaded");
    }

    const int qLen = static_cast<int>(query.size());
    const int tLen = static_cast<int>(target.size());

    std::vector<std::vector<MatrixTraceCell>> dp(
        qLen + 1,
        std::vector<MatrixTraceCell>(tLen + 1));

    dp[0][0].cost = 0.0f;
    dp[0][0].reachable = true;
    dp[0][0].prevI = -1;
    dp[0][0].prevJ = -1;
    dp[0][0].op = AlignmentOpType::Match;

    for (int i = 0; i <= qLen; ++i) {
        for (int j = 0; j <= tLen; ++j) {
            if (!dp[i][j].reachable) {
                continue;
            }

            if (i < qLen) {
                const float candidateCost = dp[i][j].cost + substitutionMatrix_.at(query[i]).at('-');
                if (BetterMatrixCandidate(candidateCost, AlignmentOpType::Deletion, dp[i + 1][j])) {
                    dp[i + 1][j].cost = candidateCost;
                    dp[i + 1][j].reachable = true;
                    dp[i + 1][j].prevI = i;
                    dp[i + 1][j].prevJ = j;
                    dp[i + 1][j].op = AlignmentOpType::Deletion;
                }
            }

            if (j < tLen) {
                const float candidateCost = dp[i][j].cost + substitutionMatrix_.at('-').at(target[j]);
                if (BetterMatrixCandidate(candidateCost, AlignmentOpType::Insertion, dp[i][j + 1])) {
                    dp[i][j + 1].cost = candidateCost;
                    dp[i][j + 1].reachable = true;
                    dp[i][j + 1].prevI = i;
                    dp[i][j + 1].prevJ = j;
                    dp[i][j + 1].op = AlignmentOpType::Insertion;
                }
            }

            if (i < qLen && j < tLen) {
                const auto op = (query[i] == target[j]) ? AlignmentOpType::Match : AlignmentOpType::Substitution;
                const float candidateCost = dp[i][j].cost + substitutionMatrix_.at(query[i]).at(target[j]);
                if (BetterMatrixCandidate(candidateCost, op, dp[i + 1][j + 1])) {
                    dp[i + 1][j + 1].cost = candidateCost;
                    dp[i + 1][j + 1].reachable = true;
                    dp[i + 1][j + 1].prevI = i;
                    dp[i + 1][j + 1].prevJ = j;
                    dp[i + 1][j + 1].op = op;
                }
            }
        }
    }

    if (!dp[qLen][tLen].reachable) {
        return std::nullopt;
    }
    if (maxCost.has_value() && dp[qLen][tLen].cost > *maxCost + kFloatEps) {
        return std::nullopt;
    }

    AlignmentResult result;
    result.distance = dp[qLen][tLen].cost;

    std::string queryAlignedRev;
    std::string targetAlignedRev;
    std::vector<AlignmentOp> opsRev;

    int i = qLen;
    int j = tLen;
    while (true) {
        const auto& cell = dp[i][j];
        if (cell.prevI < 0) {
            break;
        }

        switch (cell.op) {
            case AlignmentOpType::Match:
                queryAlignedRev.push_back(query[i - 1]);
                targetAlignedRev.push_back(target[j - 1]);
                break;
            case AlignmentOpType::Substitution:
                ++result.substitutions;
                queryAlignedRev.push_back(query[i - 1]);
                targetAlignedRev.push_back(target[j - 1]);
                opsRev.push_back({AlignmentOpType::Substitution, i, query[i - 1], target[j - 1]});
                break;
            case AlignmentOpType::Deletion:
                ++result.deletions;
                queryAlignedRev.push_back(query[i - 1]);
                targetAlignedRev.push_back('-');
                opsRev.push_back({AlignmentOpType::Deletion, i, query[i - 1], '-'});
                break;
            case AlignmentOpType::Insertion:
                ++result.insertions;
                queryAlignedRev.push_back('-');
                targetAlignedRev.push_back(target[j - 1]);
                opsRev.push_back({AlignmentOpType::Insertion, i, '-', target[j - 1]});
                break;
        }

        const int nextI = cell.prevI;
        const int nextJ = cell.prevJ;
        i = nextI;
        j = nextJ;
    }

    std::reverse(queryAlignedRev.begin(), queryAlignedRev.end());
    std::reverse(targetAlignedRev.begin(), targetAlignedRev.end());
    std::reverse(opsRev.begin(), opsRev.end());

    result.queryAligned = std::move(queryAlignedRev);
    result.targetAligned = std::move(targetAlignedRev);
    result.ops = std::move(opsRev);
    return result;
}

std::optional<Trie::AlignmentResult> Trie::AlignIndexHitWithMatrix(
        const std::string& query,
        size_t targetIndex,
        std::optional<float> maxCost) {
    if (targetIndex >= sequences_.size()) {
        return std::nullopt;
    }
    return AlignQueryToTargetWithMatrix(query, sequences_[targetIndex], maxCost);
}

std::vector<std::optional<Trie::AlignmentResult>> Trie::AlignIndexHitsWithMatrix(
        const std::string& query,
        const std::vector<std::pair<size_t, float>>& hits,
        std::optional<float> maxCost,
        std::optional<std::size_t> numThreads) {
    return RunParallelOrdered<std::optional<AlignmentResult>>(
        hits.size(),
        numThreads,
        [&](std::size_t i) {
            const auto& [targetIndex, dist] = hits[i];
            (void)dist;
            return AlignIndexHitWithMatrix(query, targetIndex, maxCost);
        });
}
