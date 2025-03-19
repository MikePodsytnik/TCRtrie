#include "Trie.h"

Trie::Trie(const vector<string>& patterns) : patterns_(patterns) {
    root_ = new TrieNode();
    BuildTrie();
}

Trie::Trie() {
    root_ = nullptr;
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
const int MAX_QUERY_LENGTH = 32;

vector<string> Trie::Search(const string& query, int maxEdits) {
    vector<string> results;
    int queryLength = query.size();
    if (queryLength > MAX_QUERY_LENGTH) {
        cerr << "Query length exceeds maximum allowed length." << endl;
        return results;
    }
    int initialRow[MAX_QUERY_LENGTH + 1];
    for (int i = 0; i <= queryLength; ++i) {
        initialRow[i] = i;
    }
    SearchRecursive(query, maxEdits, "", root_, initialRow, queryLength, results);
    return results;
}

void Trie::SearchRecursive(const string &query, int maxEdits, const string &currentPrefix, TrieNode* node, const int* prevRow, int queryLength, vector<string>& results) {
    int currentRow[MAX_QUERY_LENGTH + 1];

    memcpy(currentRow, prevRow, sizeof(int) * (queryLength + 1));
    string prefix = currentPrefix;

    if (!node->patterns_indices.empty() && currentRow[queryLength] <= maxEdits) {
        for (int patternIndex : node->patterns_indices) {
            results.push_back(patterns_[patternIndex]);
        }
    }

    int minVal = *min_element(currentRow, currentRow + queryLength + 1);
    if (minVal > maxEdits) return;

    for (int i = 0; i < node->children.size(); ++i) {
        TrieNode* child = node->children[i];
        if (child == nullptr) continue;
        char letter = 'A' + i;

        int nextRow[MAX_QUERY_LENGTH + 1];
        nextRow[0] = currentRow[0] + 1; // вставка
        for (int j = 1; j <= queryLength; ++j) {
            int cost = (query[j - 1] == letter) ? 0 : 1;
            nextRow[j] = min({ currentRow[j] + 1,       // удаление
                               nextRow[j - 1] + 1,        // вставка
                               currentRow[j - 1] + cost   // замена
                             });
        }
        SearchRecursive(query, maxEdits, prefix + letter, child, nextRow, queryLength, results);
    }
}

bool Trie::SearchAny(const string& query, int maxReplacements) {
    return SearchAnyRecursive(query, 0, root_, maxReplacements);
}


void Trie::BuildTrie() {
    for (int i = 0; i < patterns_.size(); ++i) {
        TrieNode* currentNode = root_;
        for (char c : patterns_[i]) {
            int index = c - 'A';
            if (currentNode->children[index] == nullptr) {
                currentNode->children[index] = new TrieNode();
            }
            currentNode = currentNode->children[index];
        }
        currentNode->patterns_indices.push_back(i);
    }
}

void Trie::DeleteTrie(TrieNode* node) {
    if (!node) return;
    for (TrieNode* childNode : node->children) {
        DeleteTrie(childNode);
    }
    delete node;
}

unordered_map<string, vector<string>> Trie::SearchForAll(const vector<string>& queries, int maxReplacements) {
    unordered_map<string, vector<string>> result;
    vector<future<pair<string, vector<string>>>> futures;

    for (const auto& query : queries) {
        futures.push_back(async(launch::async,
                                [this, query, maxReplacements]() -> pair<string, vector<string>> {
                                    vector<string> searchResult = this->Search(query, maxReplacements);
                                    return {query, searchResult};
                                }));
    }

    for (size_t i = 0; i < queries.size(); ++i) {
        auto searchResult = futures[i].get();
        result[searchResult.first] = searchResult.second;
    }

    return result;
}

Trie::TrieNode* Trie::CopyTrie(const TrieNode* node) {
    if (!node) return nullptr;

    TrieNode* newNode = new TrieNode();
    newNode->patterns_indices = node->patterns_indices;
    newNode->compressedEdge = node->compressedEdge;
    newNode->compressedLength = node->compressedLength;

    for (size_t i = 0; i < node->children.size(); ++i) {
        if (node->children[i]) {
            newNode->children[i] = CopyTrie(node->children[i]);
        }
    }

    return newNode;
}

bool Trie::SearchAnyRecursive(const string& query, int pos, TrieNode* currentNode, int replacementsLeft) {
    if (pos == query.size()) {
        return !currentNode->patterns_indices.empty();
    }

    char currentChar = query[pos];
    int index = currentChar - 'A';

    if (!currentNode->compressedEdge.empty()) {
        int matchLength = min((int)currentNode->compressedEdge.size(), (int)query.size() - pos);
        int mismatches = 0;

        for (int i = 0; i < matchLength; ++i) {
            if (query[pos + i] != currentNode->compressedEdge[i]) {
                ++mismatches;
                if (mismatches > replacementsLeft) {
                    return false;
                }
            }
        }

        pos += matchLength;
        if (pos == query.size()) {
            return !currentNode->patterns_indices.empty();
        }

        currentChar = query[pos];
        index = currentChar - 'A';
    }

    for (int i = 0; i < 30; ++i) {
        TrieNode* childNode = currentNode->children[i];
        if (childNode == nullptr) continue;

        if (i == index) {
            if (SearchAnyRecursive(query, pos + 1, childNode, replacementsLeft)) {
                return true;
            }
        } else if (replacementsLeft > 0) {
            if (SearchAnyRecursive(query, pos + 1, childNode, replacementsLeft - 1)) {
                return true;
            }
        }
    }

    return false;
}
