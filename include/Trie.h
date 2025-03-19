#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <queue>
#include <fstream>
#include <future>
#include <string>
#include <cstring>
#include <algorithm>

using namespace std;

class Trie {
public:
    struct TrieNode {
        vector<TrieNode*> children;
        vector<int> patterns_indices;
        string compressedEdge;
        int compressedLength = 0;

        TrieNode() : children(30, nullptr), compressedEdge(""), compressedLength(0) {}
    };

    Trie(const vector<string>& patterns);
    Trie();
    Trie(const Trie& other);
    Trie& operator=(const Trie& other);
    Trie(Trie&& other) noexcept;
    Trie& operator=(Trie&& other) noexcept;
    ~Trie();

    vector<string> Search(const string& query, int maxEdits);

    bool SearchAny(const string& query, int maxReplacements);

    unordered_map<string, vector<string>> SearchForAll(const vector<string>& queries, int maxReplacements);

private:
    static const int MAX_QUERY_LENGTH = 32;
    TrieNode* root_;
    vector<string> patterns_;

    void BuildTrie();

    void SearchRecursive(const string &query, int maxEdits, const string &currentPrefix, TrieNode* node, const int* prevRow, int queryLength, vector<string>& results);

    bool SearchAnyRecursive(const string& query, int pos, TrieNode* currentNode, int replacementsLeft);

    void DeleteTrie(TrieNode* node);

    TrieNode* CopyTrie(const TrieNode* node);
};
