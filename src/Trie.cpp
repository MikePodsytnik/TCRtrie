#include "Trie.h"
#include <cstring>

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
        nextRow[0] = currentRow[0] + 1;
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

unordered_map<string, vector<string>> Trie::SearchForAll(const vector<string>& queries, int maxReplacements) {
    unordered_map<string, vector<string>> result;
    vector<future<pair<string, vector<string>>>> futures;
    
    size_t maxConcurrent = std::thread::hardware_concurrency();
    if (maxConcurrent == 0) {
        maxConcurrent = 32;
    }

    for (const auto& query : queries) {
        if (futures.size() >= maxConcurrent) {
            auto completed = futures.front().get();
            result[completed.first] = completed.second;
            futures.erase(futures.begin());
        }

        futures.push_back(std::async(std::launch::async,
                                     [this, query, maxReplacements]() -> pair<string, vector<string>> {
                                         vector<string> searchResult = this->Search(query, maxReplacements);
                                         return {query, searchResult};
                                     }));
    }

    for (auto& f : futures) {
        auto completed = f.get();
        result[completed.first] = completed.second;
    }

    return result;
}

bool Trie::SearchAny(const string& query, int maxEdits) {
    int queryLength = query.size();
    if (queryLength > MAX_QUERY_LENGTH) {
        cerr << "Query length exceeds maximum allowed length." << endl;
        return false;
    }
    int initialRow[MAX_QUERY_LENGTH + 1];
    for (int i = 0; i <= queryLength; ++i) {
        initialRow[i] = i;
    }
    
    return SearchAnyRecursive(query, maxEdits, root_, initialRow, queryLength);
}

bool Trie::SearchAnyRecursive(const string &query, int maxEdits, TrieNode* node, const int* prevRow, int queryLength) {
    int currentRow[MAX_QUERY_LENGTH + 1];
    memcpy(currentRow, prevRow, sizeof(int) * (queryLength + 1));
    
    if (!node->patterns_indices.empty() && currentRow[queryLength] <= maxEdits) {
        return true;
    }
    
    int minVal = *min_element(currentRow, currentRow + queryLength + 1);
    if (minVal > maxEdits) return false;
    
    for (int i = 0; i < node->children.size(); ++i) {
        TrieNode* child = node->children[i];
        if (child == nullptr) continue;
        char letter = 'A' + i;
        
        int nextRow[MAX_QUERY_LENGTH + 1];
        nextRow[0] = currentRow[0] + 1;
        
        for (int j = 1; j <= queryLength; ++j) {
            int cost = (query[j - 1] == letter) ? 0 : 1;
            nextRow[j] = min({ currentRow[j] + 1,       // удаление
                               nextRow[j - 1] + 1,        // вставка
                               currentRow[j - 1] + cost   // замена/совпадение
                             });
        }
        
        if (SearchAnyRecursiveDPInPlace(query, maxEdits, child, nextRow, queryLength)) {
            return true;
        }
    }
  
    return false;
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
