#include "Trie.h"

Trie::Trie(const std::string& dataPath) {
    root_ = new TrieNode();
    LoadAIRRAndBuildTrie(dataPath);
}

Trie::Trie() : root_(nullptr) {}

Trie::Trie(Trie&& other) noexcept
        : root_(other.root_), patterns_(std::move(other.patterns_)) {
    other.root_ = nullptr;
}

Trie& Trie::operator=(Trie&& other) noexcept {
    if (this != &other) {
        DeleteTrie(root_);
        root_ = other.root_;
        patterns_ = std::move(other.patterns_);
        other.root_ = nullptr;
    }
    return *this;
}

Trie::Trie(const Trie& other) : patterns_(other.patterns_) {
    root_ = CopyTrie(other.root_);
}

Trie& Trie::operator=(const Trie& other) {
    if (this != &other) {
        DeleteTrie(root_);
        root_ = CopyTrie(other.root_);
        patterns_ = other.patterns_;
    }

    return *this;
}

Trie::~Trie() {
    DeleteTrie(root_);
}

std::vector<AIRREntity> Trie::Search(const std::string& query, int maxEdits,
                                     const std::optional<std::string>& vGeneFilter,
                                     const std::optional<std::string>& jGeneFilter) {
    std::vector<AIRREntity> results;
    int queryLength = query.size();

    if (queryLength > maxQueryLength_) {
        std::cerr << query << " :query length exceeds maximum allowed length(32)" << std::endl;
        return results;
    }

    int initialRow[maxQueryLength_ + 1];
    for (int i = 0; i <= queryLength; ++i) {
        initialRow[i] = i;
    }

    SearchRecursive(query, maxEdits, "", root_, initialRow, queryLength, results, vGeneFilter, jGeneFilter);
    return results;
}

void Trie::SearchRecursive(const std::string &query, int maxEdits, const std::string &currentPrefix,
                           TrieNode* node, const int* prevRow, int queryLength,
                           std::vector<AIRREntity>& results,
                           const std::optional<std::string>& vGeneFilter,
                           const std::optional<std::string>& jGeneFilter) {
    int currentRow[maxQueryLength_ + 1];
    memcpy(currentRow, prevRow, sizeof(int) * (queryLength + 1));

    if (!node->patternsIndices.empty() && currentRow[queryLength] <= maxEdits) {
        for (int patternIndex : node->patternsIndices) {
            bool vMatch = !vGeneFilter || patterns_[patternIndex].vGene == *vGeneFilter;
            bool jMatch = !jGeneFilter || patterns_[patternIndex].jGene == *jGeneFilter;
            if (vMatch && jMatch) {
                results.push_back(patterns_[patternIndex]);
            }
        }
    }

    int minVal = *std::min_element(currentRow, currentRow + queryLength + 1);
    if (minVal > maxEdits) return;

    for (int i = 0; i < node->children.size(); ++i) {
        TrieNode* child = node->children[i];
        if (child == nullptr) continue;
        char letter = 'A' + i;

        int nextRow[maxQueryLength_ + 1];
        nextRow[0] = currentRow[0] + 1;
        for (int j = 1; j <= queryLength; ++j) {
            int cost = (query[j - 1] == letter) ? 0 : 1;
            nextRow[j] = std::min({ currentRow[j] + 1,
                                    nextRow[j - 1] + 1,
                                    currentRow[j - 1] + cost
                                  });
        }
        SearchRecursive(query, maxEdits, currentPrefix + letter, child, nextRow, queryLength, results, vGeneFilter, jGeneFilter);
    }
}

std::vector<AIRREntity> Trie::SearchWithScore(const std::string& query, int maxScore,
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

    int initialRow[maxQueryLength_ + 1];
    for (int i = 0; i <= queryLength; ++i) {
        initialRow[i] = i * depletionCost;
    }

    SearchRecursiveScore(query, maxScore, "", root_, initialRow, queryLength,
                         results, vGeneFilter, jGeneFilter);

    return results;
}

void Trie::SearchRecursiveScore(const std::string &query, int maxScore, const std::string &currentPrefix,
                                TrieNode* node, const int* prevRow, int queryLength,
                                std::vector<AIRREntity>& results,
                                const std::optional<std::string>& vGeneFilter,
                                const std::optional<std::string>& jGeneFilter) {
    int currentRow[maxQueryLength_ + 1];
    memcpy(currentRow, prevRow, sizeof(int) * (queryLength + 1));

    if (!node->patternsIndices.empty() && currentRow[queryLength] <= maxScore) {
        for (int patternIndex : node->patternsIndices) {
            bool vMatch = !vGeneFilter || patterns_[patternIndex].vGene == *vGeneFilter;
            bool jMatch = !jGeneFilter || patterns_[patternIndex].jGene == *jGeneFilter;

            if (vMatch && jMatch) {
                results.push_back(patterns_[patternIndex]);
            }
        }
    }

    int minVal = *std::min_element(currentRow, currentRow + queryLength + 1);
    if (minVal > maxScore) return;

    for (int i = 0; i < node->children.size(); ++i) {
        TrieNode* child = node->children[i];
        if (!child) continue;
        char letter = 'A' + i;

        int nextRow[maxQueryLength_ + 1];
        nextRow[0] = currentRow[0] + insertionCost;

        for (int j = 1; j <= queryLength; ++j) {
            char queryChar = query[j - 1];
            int subCost = substitutionMatrix_.at(queryChar).at(letter);

            nextRow[j] = std::min({
                                          currentRow[j] + insertionCost,
                                          nextRow[j - 1] + depletionCost,
                                          currentRow[j - 1] + subCost
                                  });
        }

        SearchRecursiveScore(query, maxScore, currentPrefix + letter, child, nextRow, queryLength,
                             results, vGeneFilter, jGeneFilter);
    }
}

std::unordered_map<std::string, std::vector<AIRREntity>> Trie::SearchForAll(
        const std::vector<std::string>& queries,
        int maxEdits,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter) {

    std::unordered_map<std::string, std::vector<AIRREntity>> result;
    std::vector<std::future<std::pair<std::string, std::vector<AIRREntity>>>> futures;

    std::size_t maxConcurrent = 10 * std::thread::hardware_concurrency();

    for (std::size_t i = 0; i < queries.size(); ++i) {
        const std::string& query = queries[i];

        futures.emplace_back(std::async(std::launch::async,
                                        [this, query, maxEdits]() -> std::pair<std::string, std::vector<AIRREntity>> {
                                            return { query, this->Search(query, maxEdits) };
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

std::unordered_map<std::string, std::vector<AIRREntity>> Trie::SearchForAllWithScore(
        const std::vector<std::string>& queries,
        int maxScore,
        const std::optional<std::string>& vGeneFilter,
        const std::optional<std::string>& jGeneFilter) {

    std::unordered_map<std::string, std::vector<AIRREntity>> result;
    std::vector<std::future<std::pair<std::string, std::vector<AIRREntity>>>> futures;

    std::size_t maxConcurrent = 10 * std::thread::hardware_concurrency();

    for (std::size_t i = 0; i < queries.size(); ++i) {
        const std::string& query = queries[i];

        futures.push_back(std::async(std::launch::async,
                                     [this, query, maxScore]() {
                                         return std::make_pair(query, this->SearchWithScore(query, maxScore));
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
    int initialRow[maxQueryLength_ + 1];
    for (int i = 0; i <= queryLength; ++i) {
        initialRow[i] = i;
    }

    return SearchAnyRecursive(query, maxEdits, root_, initialRow, queryLength);
}

bool Trie::SearchAnyRecursive(const std::string &query, int maxEdits,
                              TrieNode* node, const int* prevRow, int queryLength) {
    int currentRow[maxQueryLength_ + 1];
    memcpy(currentRow, prevRow, sizeof(int) * (queryLength + 1));

    if (!node->patternsIndices.empty() && currentRow[queryLength] <= maxEdits) {
        for (int patternIndex : node->patternsIndices) {
            return true;
        }
    }

    int minVal = *std::min_element(currentRow, currentRow + queryLength + 1);
    if (minVal > maxEdits) return false;

    for (int i = 0; i < node->children.size(); ++i) {
        TrieNode* child = node->children[i];
        if (child == nullptr) continue;
        char letter = 'A' + i;

        int nextRow[maxQueryLength_ + 1];
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

void Trie::LoadAIRRAndBuildTrie(const std::string& dataPath) {
    auto entries = ParseAIRR(dataPath);
    int patternIndex = 0;

    for (const auto& entry : entries) {
        patterns_.push_back(entry);

        TrieNode* node = root_;
        for (char c : entry.junctionAA) {
            if (c < 'A' || c > 'Z') continue;
            int idx = c - 'A';
            if (!node->children[idx]) node->children[idx] = new TrieNode();
            node = node->children[idx];
        }
        node->patternsIndices.push_back(patternIndex++);
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
    newNode->patternsIndices = node->patternsIndices;

    for (size_t i = 0; i < node->children.size(); ++i) {
        if (node->children[i]) {
            newNode->children[i] = CopyTrie(node->children[i]);
        }
    }

    return newNode;
}

void Trie::LoadSubstitutionMatrix(const std::string& matrixPath) {
    std::unordered_map<char, std::unordered_map<char, float>> substitutionMatrix;

    std::ifstream file(matrixPath);
    if (!file) {
        std::cerr << "Error: Unable to open substitution matrix file." << std::endl;
        return;
    }

    std::vector<char> letters;
    std::string line;
    getline(file, line);
    std::istringstream headerStream(line);
    char letter;
    while (headerStream >> letter) {
        letters.push_back(letter);
    }

    std::unordered_map<char, std::unordered_map<char, float>> rawScores;
    float maxScore = -1e9;

    for (char rowLetter : letters) {
        std::unordered_map<char, float> rowMap;
        file >> letter;

        for (char colLetter : letters) {
            float score;
            file >> score;
            rowMap[colLetter] = score;
            if (score > maxScore) {
                maxScore = score;
            }
        }

        rawScores[rowLetter] = rowMap;
    }

    for (const auto& rowPair : rawScores) {
        char rowLetter = rowPair.first;
        const auto& rowMap = rowPair.second;
        std::unordered_map<char, float> costRow;

        for (const auto& colPair : rowMap) {
            char colLetter = colPair.first;
            float score = colPair.second;
            float cost = maxScore - score;
            costRow[colLetter] = cost;
        }

        substitutionMatrix[rowLetter] = costRow;
    }

    substitutionMatrix_ = substitutionMatrix;
    useSubstitutionMatrix_ = true;
}

void Trie::SetMaxQueryLength(int newMaxQueryLength) {
    maxQueryLength_ = newMaxQueryLength;
}

bool Trie::IsStandardAminoAcidSequence(const std::string& seq) const {
    static const std::unordered_set<char> validAminoAcids = {
            'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I',
            'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
    };

    for (char c : seq) {
        if (!validAminoAcids.count(c)) {
            return false;
        }
    }
    return true;
}