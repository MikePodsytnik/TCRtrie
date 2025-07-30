#include "Trie.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <sstream>

Trie::Trie(const std::string& dataPath) {
    root_ = new TrieNode();
    LoadAIRR(dataPath);
    BuildTrie();
}

Trie::Trie(const std::vector<std::string>& sequences)
        : root_(new TrieNode()), sequences_(sequences)
{
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
          jGenes_(other.jGenes_)
{
    root_ = CopyTrie(other.root_);
}

Trie::Trie(Trie&& other) noexcept
        : root_(other.root_),
          maxQueryLength_(other.maxQueryLength_),
          useSubstitutionMatrix_(other.useSubstitutionMatrix_),
          substitutionMatrix_(std::move(other.substitutionMatrix_)),
          sequences_(std::move(other.sequences_)),
          vGenes_(std::move(other.vGenes_)),
          jGenes_(std::move(other.jGenes_))
{
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

std::vector<Trie::Stat> Trie::PruneStats(const std::vector<Trie::Stat>& stats) {
    std::vector<Trie::Stat> res;
    for (auto const& st : stats) {
        bool dominated = false;
        for (auto it = res.begin(); it != res.end(); ) {
            auto const& ex = *it;
            if (ex.distance <= st.distance
                && ex.insertion <= st.insertion
                && ex.deletion <= st.deletion
                && ex.substitution <= st.substitution) {
                dominated = true;
                break;
            }
            if (st.distance <= ex.distance
                && st.insertion <= ex.insertion
                && st.deletion <= ex.deletion
                && st.substitution <= ex.substitution) {
                it = res.erase(it);
            } else {
                ++it;
            }
        }
        if (!dominated)
            res.push_back(st);
    }
    return res;
}

std::vector<Trie::Stat> Trie::DetailedLevenshteinAll(
        const std::string& s,
        const std::string& t,
        int maxEdits
) {
    int n = s.size(), m = t.size();
    std::vector<std::vector<std::vector<Trie::Stat>>> dp(
            n+1, std::vector<std::vector<Trie::Stat>>(m+1));

    dp[0][0].push_back({0,0,0,0});
    for (int i = 1; i <= n; ++i) {
        dp[i][0].push_back({i, 0, i, 0});
    }
    for (int j = 1; j <= m; ++j) {
        dp[0][j].push_back({j, j, 0, 0});
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            std::vector<Trie::Stat> cand;
            for (auto const& st : dp[i-1][j]) {
                Stat ns = st;
                ns.distance++; ns.deletion++;
                if (ns.distance <= maxEdits) cand.push_back(ns);
            }
            for (auto const& st : dp[i][j-1]) {
                Stat ns = st;
                ns.distance++; ns.insertion++;
                if (ns.distance <= maxEdits) cand.push_back(ns);
            }
            int cost = (s[i-1] == t[j-1] ? 0 : 1);
            for (auto const& st : dp[i-1][j-1]) {
                Stat ns = st;
                ns.distance += cost;
                ns.substitution += cost;
                if (ns.distance <= maxEdits) cand.push_back(ns);
            }
            dp[i][j] = PruneStats(cand);
        }
    }
    return dp[n][m];
}

std::vector<AIRREntity> Trie::SearchAIRR(const std::string& query,
                                         int maxSubstitution,
                                         int maxInsertion,
                                         int maxDeletion,
                                         const std::optional<std::string>& vGeneFilter,
                                         const std::optional<std::string>& jGeneFilter) {
    int maxEdits = maxSubstitution + maxInsertion + maxDeletion;
    std::vector<AIRREntity> results;
    int queryLength = query.size();

    if (queryLength > maxQueryLength_) {
        std::cerr << query << " :query length exceeds maximum allowed length(" << maxQueryLength_ << ")" << std::endl;
        return results;
    }

    std::vector<int> initialRow(maxQueryLength_ + 1);
    for (int i = 0; i <= queryLength; ++i) {
        initialRow[i] = i;
    }

    SearchRecursiveAIRR(query, maxEdits, root_, initialRow, queryLength, results, vGeneFilter, jGeneFilter);
    std::vector<AIRREntity> finalResult;
    for (const auto& candidate : results) {
        auto allStats = DetailedLevenshteinAll(query, candidate.junctionAA, maxEdits);
        bool ok = false;
        for (auto& st : allStats) {
            if (st.substitution <= maxSubstitution
                && st.insertion <= maxInsertion
                && st.deletion <= maxDeletion) {
                ok = true;
                break;
            }
        }
        if (ok) {
            finalResult.push_back(candidate);
        }
    }

    return finalResult;
}

void Trie::SearchRecursiveAIRR(const std::string& query, int maxEdits,
                               TrieNode* node, std::vector<int>& prevRow, int queryLength,
                               std::vector<AIRREntity>& results,
                               const std::optional<std::string>& vGeneFilter,
                               const std::optional<std::string>& jGeneFilter) {
    std::vector<int> currentRow(maxQueryLength_ + 1);
    std::copy(prevRow.begin(), prevRow.begin() + queryLength + 1, currentRow.begin());

    if (!node->indices.empty() && currentRow[queryLength] <= maxEdits) {
        for (int index : node->indices) {
            bool vMatch = !vGeneFilter || vGenes_[index] == *vGeneFilter;
            bool jMatch = !jGeneFilter || jGenes_[index] == *jGeneFilter;
            if (vMatch && jMatch) {
                results.emplace_back(sequences_[index],
                                     vGenes_[index],
                                     jGenes_[index],
                                     currentRow[queryLength]);
            }
        }
    }

    int minVal = *std::min_element(currentRow.begin(), currentRow.begin() + queryLength + 1);
    if (minVal > maxEdits) return;

    for (int i = 0; i < node->children.size(); ++i) {
        TrieNode* child = node->children[i];
        if (child == nullptr) continue;
        char letter = 'A' + i;

        std::vector<int> nextRow(maxQueryLength_ + 1);
        nextRow[0] = currentRow[0] + 1;
        for (int j = 1; j <= queryLength; ++j) {
            int cost = (query[j - 1] == letter) ? 0 : 1;
            nextRow[j] = std::min({ currentRow[j] + 1,
                                    nextRow[j - 1] + 1,
                                    currentRow[j - 1] + cost
                                  });
        }
        SearchRecursiveAIRR(query, maxEdits, child, nextRow, queryLength, results, vGeneFilter, jGeneFilter);
    }
}

std::vector<AIRREntity> Trie::SearchWithMatrix(const std::string& query, float maxCost,
                                               const std::optional<std::string>& vGeneFilter,
                                               const std::optional<std::string>& jGeneFilter) {
    std::vector<AIRREntity> results;
    int queryLength = query.size();

    if (!useSubstitutionMatrix_) {
        std::cerr << "No substitution matrix is entered, only Levenshtein distance search is available" << std::endl;
        return results;
    }

    if (queryLength > maxQueryLength_) {
        std::cerr << "Query length exceeds maximum allowed length." << std::endl;
        return results;
    }

    std::vector<float> initialRow(maxQueryLength_ + 1);
    initialRow[0] = 0;
    for (int i = 1; i <= queryLength; ++i) {
        initialRow[i] = initialRow[i-1] + substitutionMatrix_.at('-').at(query[i-1]);
    }
    SearchRecursiveCost(query, maxCost, root_, initialRow, queryLength,
                        results, vGeneFilter, jGeneFilter);

    return results;
}

void Trie::SearchRecursiveCost(const std::string& query, float maxCost,
                               TrieNode* node, std::vector<float>& prevRow, int queryLength,
                               std::vector<AIRREntity>& results,
                               const std::optional<std::string>& vGeneFilter,
                               const std::optional<std::string>& jGeneFilter) {
    std::vector<float> currentRow(maxQueryLength_ + 1);
    std::copy(prevRow.begin(), prevRow.begin() + queryLength + 1, currentRow.begin());


    if (!node->indices.empty() && (currentRow[queryLength] <= maxCost)) {
        for (int index : node->indices) {
            bool vMatch = !vGeneFilter || vGenes_[index] == *vGeneFilter;
            bool jMatch = !jGeneFilter || jGenes_[index] == *jGeneFilter;
            if (vMatch && jMatch) {
                results.emplace_back(sequences_[index],
                                     vGenes_[index],
                                     jGenes_[index],
                                     currentRow[queryLength]);
            }
        }
    }

    for (int i = 0; i < node->children.size(); ++i) {
        TrieNode* child = node->children[i];
        if (!child) continue;
        char letter = 'A' + i;

        std::vector<float> nextRow(maxQueryLength_ + 1);
        nextRow[0] = currentRow[0] + substitutionMatrix_.at('-').at(letter);
        float minVal = nextRow[0];

        float depletionCost = substitutionMatrix_.at('-').at(letter);
        for (int j = 1; j <= queryLength; ++j) {
            char queryChar = query[j - 1];
            float subCost = substitutionMatrix_.at(queryChar).at(letter);
            float insertionCost = substitutionMatrix_.at('-').at(queryChar);

            nextRow[j] = std::min({
                                          currentRow[j] + depletionCost,
                                          nextRow[j - 1] + insertionCost,
                                          currentRow[j - 1] + subCost
                                  });
            minVal = std::min(minVal, nextRow[j]);
        }

        if (minVal > maxCost) continue;

        SearchRecursiveCost(query, maxCost, child, nextRow, queryLength,
                            results, vGeneFilter, jGeneFilter);
    }
}

std::vector<std::string> Trie::Search(const std::string& query, int maxEdits) {
    std::vector<std::string> results;
    int queryLength = query.size();
    if (queryLength > maxQueryLength_) {
        std::cerr << "Query length exceeds maximum allowed length." << std::endl;
        return results;
    }
    std::vector<int> initialRow(maxQueryLength_ + 1);
    for (int i = 0; i <= queryLength; ++i) {
        initialRow[i] = i;
    }
    SearchRecursive(query, maxEdits, "", root_, initialRow, queryLength, results);

    return results;
}

void Trie::SearchRecursive(const std::string& query, int maxEdits, const std::string& currentPrefix,
                           TrieNode* node, std::vector<int>& prevRow, int queryLength, std::vector<std::string>& results) {
    std::vector<int> currentRow(maxQueryLength_ + 1);

    std::copy(prevRow.begin(), prevRow.begin() + queryLength + 1, currentRow.begin());
    std::string prefix = currentPrefix;

    if (!node->indices.empty() && currentRow[queryLength] <= maxEdits) {
        for (int index : node->indices) {
            results.push_back(sequences_[index]);
        }
    }

    int minVal = *std::min_element(currentRow.begin(), currentRow.begin() + queryLength + 1);
    if (minVal > maxEdits) return;

    for (int i = 0; i < node->children.size(); ++i) {
        TrieNode* child = node->children[i];
        if (child == nullptr) continue;
        char letter = 'A' + i;

        std::vector<int> nextRow(maxQueryLength_ + 1);
        nextRow[0] = currentRow[0] + 1;
        for (int j = 1; j <= queryLength; ++j) {
            int cost = (query[j - 1] == letter) ? 0 : 1;
            nextRow[j] = std::min({ currentRow[j] + 1,
                                    nextRow[j - 1] + 1,
                                    currentRow[j - 1] + cost
                                  });
        }
        SearchRecursive(query, maxEdits, prefix + letter, child, nextRow, queryLength, results);
    }
}

std::unordered_map<std::string, std::vector<std::string>> Trie::Search(const std::vector<std::string>& queries,
                                                                       int maxEdits) {
    std::unordered_map<std::string, std::vector<std::string>> result;
    std::vector<std::future<std::pair<std::string, std::vector<std::string>>>> futures;

    std::size_t maxConcurrent = 10 * std::thread::hardware_concurrency();

    for (std::size_t i = 0; i < queries.size(); ++i) {
        const std::string& query = queries[i];

        futures.emplace_back(std::async(std::launch::async,
                                        [this, query, maxEdits]() -> std::pair<std::string, std::vector<std::string>> {
                                            return { query, this->Search(query, maxEdits) };
                                        }));

        if (futures.size() >= maxConcurrent || i == queries.size() - 1) {
            for (std::future<std::pair<std::string, std::vector<std::string>>>& fut : futures) {
                std::pair<std::string, std::vector<std::string>> completed = fut.get();
                result[std::move(completed.first)] = std::move(completed.second);
            }
            futures.clear();
        }
    }

    return result;
}

std::unordered_map<std::string, std::vector<AIRREntity>> Trie::SearchForAll(
        const std::vector<std::string>& queries,
        int maxSubstitution,
        int maxInsertion,
        int maxDeletion,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter) {

    std::unordered_map<std::string, std::vector<AIRREntity>> result;
    std::vector<std::future<std::pair<std::string, std::vector<AIRREntity>>>> futures;

    std::size_t maxConcurrent = 10 * std::thread::hardware_concurrency();

    for (std::size_t i = 0; i < queries.size(); ++i) {
        const std::string& query = queries[i];

        futures.emplace_back(std::async(std::launch::async,
                                        [this,
                                                query,
                                                maxSubstitution,
                                                maxInsertion,
                                                maxDeletion,
                                                vGeneFilter,
                                                jGeneFilter]() -> std::pair<std::string, std::vector<AIRREntity>> {
                                            return { query,
                                                     this->SearchAIRR(query,
                                                                      maxSubstitution,
                                                                      maxInsertion,
                                                                      maxDeletion,
                                                                      vGeneFilter,
                                                                      jGeneFilter) };
                                        }));

        if (futures.size() >= maxConcurrent || i == queries.size() - 1) {
            for (std::future<std::pair<std::string, std::vector<AIRREntity>>>& fut : futures) {
                std::pair<std::string, std::vector<AIRREntity>> completed = fut.get();
                result[std::move(completed.first)] = std::move(completed.second);
            }
            futures.clear();
        }
    }

    return result;
}

std::unordered_map<std::string, std::vector<AIRREntity>> Trie::SearchForAllWithMatrix(
        const std::vector<std::string>& queries,
        float maxCost,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter) {

    std::unordered_map<std::string, std::vector<AIRREntity>> result;
    std::vector<std::future<std::pair<std::string, std::vector<AIRREntity>>>> futures;

    std::size_t maxConcurrent = 10 * std::thread::hardware_concurrency();

    for (std::size_t i = 0; i < queries.size(); ++i) {
        const std::string& query = queries[i];

        futures.push_back(std::async(std::launch::async,
                                     [this, query, maxCost, vGeneFilter, jGeneFilter]() {
                                         return std::make_pair(query, this->SearchWithMatrix(query,
                                                                                             maxCost,
                                                                                             vGeneFilter,
                                                                                             jGeneFilter));
                                     }));

        if (futures.size() >= maxConcurrent || i == queries.size() - 1) {
            for (std::future<std::pair<std::string, std::vector<AIRREntity>>>& fut : futures) {
                std::pair<std::string, std::vector<AIRREntity>> completed = fut.get();
                result[std::move(completed.first)] = std::move(completed.second);
            }
            futures.clear();
        }
    }

    for (auto& fut : futures) {
        auto [query, matches] = fut.get();
        result[query] = std::move(matches);
    }

    return result;
}

bool Trie::SearchAny(const std::string& query, int maxEdits) {
    int queryLength = query.size();
    if (queryLength > maxQueryLength_) {
        std::cerr << "Query length exceeds maximum allowed length." << std::endl;
        return false;
    }
    std::vector<int> initialRow(maxQueryLength_ + 1);
    for (int i = 0; i <= queryLength; ++i) {
        initialRow[i] = i;
    }

    return SearchAnyRecursive(query, maxEdits, root_, initialRow, queryLength);
}

bool Trie::SearchAnyRecursive(const std::string& query, int maxEdits,
                              TrieNode* node, std::vector<int>& prevRow, int queryLength) {
    std::vector<int> currentRow(maxQueryLength_ + 1);
    std::copy(prevRow.begin(), prevRow.begin() + queryLength + 1, currentRow.begin());

    if (!node->indices.empty() && currentRow[queryLength] <= maxEdits) {
        for (int index : node->indices) {
            return true;
        }
    }

    int minVal = *std::min_element(currentRow.begin(), currentRow.begin() + queryLength + 1);
    if (minVal > maxEdits) return false;

    for (int i = 0; i < node->children.size(); ++i) {
        TrieNode* child = node->children[i];
        if (child == nullptr) continue;
        char letter = 'A' + i;

        std::vector<int> nextRow(maxQueryLength_ + 1);
        nextRow[0] = currentRow[0] + 1;

        for (int j = 1; j <= queryLength; ++j) {
            int cost = (query[j - 1] == letter) ? 0 : 1;
            nextRow[j] = std::min({
                                          currentRow[j] + 1,
                                          nextRow[j - 1] + 1,
                                          currentRow[j - 1] + cost
                                  });
        }

        if (SearchAnyRecursive(query, maxEdits, child, nextRow, queryLength)) {
            return true;
        }
    }

    return false;
}

void Trie::LoadAIRR(const std::string& dataPath) {
    auto entries = ParseAIRR(dataPath);
    for (auto& e : entries) {
        sequences_.push_back(e.junctionAA);
        vGenes_.push_back(e.vGene);
        jGenes_.push_back(e.jGene);
    }
}

void Trie::BuildTrie() {
    for (int idx = 0; idx < sequences_.size(); ++idx) {
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
    for (TrieNode* childNode : node->children) {
        DeleteTrie(childNode);
    }
    delete node;
}

Trie::TrieNode* Trie::CopyTrie(const TrieNode* node) {
    if (!node) return nullptr;

    TrieNode* newNode = new TrieNode();
    newNode->indices = node->indices;

    for (size_t i = 0; i < node->children.size(); ++i) {
        if (node->children[i]) {
            newNode->children[i] = CopyTrie(node->children[i]);
        }
    }

    return newNode;
}

void Trie::LoadSubstitutionMatrix(const std::string& matrixPath) {
    std::ifstream file(matrixPath);
    if (!file) {
        std::cerr << "Cannot open matrix\n";
        return;
    }

    std::vector<char> letters;
    std::string line;
    std::getline(file, line);
    std::istringstream hs(line);
    char letter;
    while (hs >> letter) {
        letters.push_back(letter);
    }

    std::unordered_map<char, std::unordered_map<char, float>> rawScores;
    bool isCostMatrix = true;
    rawScores['-']['-'] = fabs(deletionScore_);
    for (char r : letters) {
        rawScores[r]['-']=deletionScore_;
        rawScores['-'][r]=deletionScore_;
        file >> letter;
        for (char c : letters) {
            float v;
            file >> v;
            rawScores[r][c] = v;
            rawScores[c][r] = v;
            if (v > 1e-6f) { isCostMatrix = false; }
        }
    }
    deletionScore_ = rawScores['-']['-'];

    letters.push_back('-');

    substitutionMatrix_.clear();

    if (isCostMatrix) {
        substitutionMatrix_ = rawScores;
    } else {
        for (char r : letters) {
            for (char c : letters) {
                float cost = (rawScores[r][r] + rawScores[c][c]) * 0.5f
                             - rawScores[r][c];
                substitutionMatrix_[r][c] = cost;
            }
        }
    }

    useSubstitutionMatrix_ = true;

    std::cout << "Substitution-Score Matrix:" << std::endl;
    PrintMatrix();
}

void Trie::PrintMatrix() {
    std::vector<char> keys;
    keys.reserve(substitutionMatrix_.size());
    for (const auto& kv : substitutionMatrix_) {
        keys.push_back(kv.first);
    }
    std::sort(keys.begin(), keys.end());

    std::cout << std::setw(4) << "";
    for (char c : keys) {
        std::cout << std::setw(6) << c;
    }
    std::cout << "\n";

    for (char row : keys) {
        std::cout << std::setw(4) << row;
        for (char col : keys) {
            float value = substitutionMatrix_[row][col];
            std::cout
                    << std::setw(6)
                    << std::fixed << std::setprecision(2)
                    << value;
        }
        std::cout << "\n";
    }
}

void Trie::SetMaxQueryLength(int newMaxQueryLength) {
    maxQueryLength_ = newMaxQueryLength;
}

void Trie::UpdateSubstitutionMatrix(float deletionScore) {
    std::vector<char> keys;
    keys.reserve(substitutionMatrix_.size());
    for (const auto& kv : substitutionMatrix_) {
        keys.push_back(kv.first);
    }
    for (auto c : keys) {
        if (c == '-') continue;

        substitutionMatrix_[c]['-'] -= deletionScore_ * 0.5f;
        substitutionMatrix_['-'][c] -= deletionScore_ * 0.5f;
        substitutionMatrix_[c]['-'] += fabs(deletionScore) * 0.5f;
        substitutionMatrix_['-'][c] += fabs(deletionScore) * 0.5f;
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