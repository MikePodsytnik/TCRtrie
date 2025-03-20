# TCRtrie: Trie-Based Neighboring CDR3 Sequences Search

TCRtrie is a C++ library and command-line tool for approximate matching of T-cell receptor (TCR) sequences using a trie (prefix tree) data structure. This tool is designed for efficient similarity search in large TCR repertoire datasets, enabling researchers to quickly identify sequences that match a given query within a specified error threshold.

## Overview

Defining and measuring similarity between TCR sequences is challenging due to the immense diversity of CDR3 sequences. TCRtrie addresses this problem by:
- Building a **trie** from a set of CDR3 sequences (patterns).
- Employing a dynamic programming approach (based on **Levenshtein** distance) to perform approximate searches.
- Utilizing multithreading to accelerate search operations across large datasets.
- Support for substitution scoring matrices (e.g. BLOSUM) for sequences at a certain score distance.

## Functional

### Search
**Signature:** `vector<string> Search(const string& query, int maxEdits)`

**Description:** Searches a Trie and returns a list of sequences whose Levenshtein distance does not exceed `maxEdits`.

### SearchAny
**Signature:** `bool SearchAny(const string& query, int maxEdits)`

**Description:** Returns `true` if at least one sequence satisfies the approximate match condition with the given query, otherwise — `false`.

### SearchForAll
**Signature:** `unordered_map<string, vector<string>> SearchForAll(const vector<string>& queries, int maxEdits)`

**Description:** Performs a parallel search for each query in the list and returns a mapping of the query to a list of matching sequences.

## How It Works

1. **Trie Construction:**  
   The trie is built from a list of TCR sequences (patterns). Each node in the trie contains:
    - A fixed-size array of child pointers (for letters 'A' to 'Z', adjusted for the alphabet used).
    - A list of indices corresponding to the patterns that terminate at that node.

2. **Approximate Search:**  
   When a query is executed:
    - The algorithm initializes a row of edit distances.
    - A recursive search traverses the trie, computing the Levenshtein distance for each node. As a result, instead of the conventional two-dimensional dynamic programming matrix, a branched, multi-dimensional matrix is obtained.
    - If a node’s computed distance is within the allowed maximum edits, the corresponding CDR3 sequences are returned.

3. **Multithreaded Batch Processing:**  
   The `SearchForAll` function uses C++ standard threading (`std::async` and `std::future`) to process multiple queries in parallel, improving performance on multi-core systems.

## Installation

### Requirements
- **C++ Compiler:** Must support at least C++17.
- **CMake:** For building the project.

### Building the Project

Clone the repository and build the project using CMake:

```sh
git clone https://github.com/yourusername/TCRtrie.git
cd TCRtrie
mkdir build && cd build
cmake ..
make
```

## Usage

**The provided main.cpp demonstrates how to:**

- Load a list of CDR3 sequences from a file (e.g., VDJdb.txt).
- Build the trie from these sequences.
- Execute searches (including multithreaded batch search) with a maximum allowed edit distance.

***Note:*** In the example provided in main.cpp, the program loads approximately 92,000 sequences from VDJdb and performs a self-search with a maximum edit distance of 3. On a MacBook Air (M1), this operation completed in approximately 21 minutes and 47 seconds. Actual performance will vary based on hardware and dataset size.

## Input Data Format
The input file should be a plain text file where each line represents one CDR3 sequence. Please ensure that:

- Sequences are written in uppercase letters (A-Z) representing amino acids.
- The file does not contain empty lines or extraneous characters.

Example (`VDJdb.txt`):
```
CASSLGTDTQYF
CASSIRSSYEQYF
CASSWGGGSHYGYTF
...
```

## Contributing
If you encounter any bugs or have suggestions for improvements, please create an issue or submit a pull request on GitHub.
