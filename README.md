# TCRtrie: Trie-Based Neighboring CDR3 Sequences Search

TCRtrie is a C++ library and command-line tool for approximate matching of T-cell receptor (TCR) sequences using a trie (prefix tree) data structure. This tool is designed for efficient similarity search in large TCR repertoire datasets, enabling researchers to quickly identify sequences that match a given query within a specified error threshold.

## Overview

Defining and measuring similarity between TCR sequences is challenging due to the immense diversity of CDR3 sequences. TCRtrie addresses this problem by:
- Building a **trie** from a set of CDR3 sequences (AIRR format).
- Employing a dynamic programming approach (based on **Levenshtein** distance) or **substitution matrix** scoring to perform approximate searches.
- Utilizing multithreading to accelerate batch search operations.
- Supporting gene-based filtering (e.g., `v_call`, `j_call`) via CLI.

## Functional

### Search
**Signature:** `vector<AIRREntity> Search(const string& query, int maxEdits, optional<string> vGene = {}, optional<string> jGene = {})`

**Description:** Searches a Trie and returns a list of matching sequences whose Levenshtein distance does not exceed `maxEdits`. Can optionally filter by V and J genes.

### SearchWithScore
**Signature:** `vector<AIRREntity> SearchWithScore(const string& query, int maxScore, optional<string> vGene = {}, optional<string> jGene = {})`

**Description:** Searches using a substitution matrix with score threshold `maxScore`. Optional filtering by V and J genes.

### SearchAny
**Signature:** `bool SearchAny(const string& query, int maxEdits)`

**Description:** Returns `true` if at least one sequence satisfies the approximate match condition with the given query.

### SearchForAll
**Signature:** `unordered_map<string, vector<AIRREntity>> SearchForAll(const vector<string>& queries, int maxEdits, optional<string> vGene = {}, optional<string> jGene = {})`

**Description:** Performs multithreaded search for all queries with Levenshtein distance. Optional gene filtering.

### SearchForAllWithScore
**Signature:** `unordered_map<string, vector<AIRREntity>> SearchForAllWithScore(const vector<string>& queries, int maxScore, optional<string> vGene = {}, optional<string> jGene = {})`

**Description:** Performs multithreaded search using a substitution matrix and score threshold. Filters supported.

### LoadSubstitutionMatrix
**Signature:** `void LoadSubstitutionMatrix(const string& matrixPath)`

**Description:** Loads a substitution matrix and converts it to a cost matrix for use in score-based search.

### SetMaxQueryLength
**Signature:** `void SetMaxQueryLength(int newMaxQueryLength)`

**Description:** Sets the maximum allowed query length (default is 32).
### CLI Interface

The project includes a command-line tool built with [CLI11](https://github.com/CLIUtils/CLI11). Example usage:

```sh
./TCRtrie \
  -i data/vdjdb.tsv \
  -q CASSLGTDGYTF \
  --n_edits 3 \
  --v_gene TRBV7-9 \
  --j_gene TRBJ2-1 \
  -o output/
```

| Flag | Description                                                                  |
|------|------------------------------------------------------------------------------|
| `-i, --input <path>` | Path to the AIRR TSV file containing the repertoire to search (**required**) |
| `-q, --query <sequence>` | Single query sequence                                                        |
| `--input_queries <path>` | AIRR TSV file with multiple queries (batch search)                           |
| `--n_edits <int>` | Max allowed Levenshtein distance                                             |
| `--matrix_search <path>` | Path to substitution matrix file                                             |
| `--score_radius <float>` | Score threshold when using matrix search                                     |
| `--v_gene <name>` | Optional filter by V-gene name                                               |
| `--j_gene <name>` | Optional filter by J-gene name                                               |
| `-o, --output <dir>` | Output folder (default: current directory)                                   |

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
### Input Format

Input files must conform to the AIRR standard (TSV) and contain at least the column `junction_aa`. Columns `v_call` and `j_call` are optional, but if any line includes one of them, all lines must include it.

Example:
```
junction_aa	v_call	j_call
CASSLGTDGYTF	TRBV7-9	TRBJ2-1
...
```

### Substitution Matrix Format

If you are using `--matrix_search`, the matrix file should follow a **BLOSUM-like format**:
- The first line contains amino acid letters (space-separated).
- Each following line begins with a row letter and contains float substitution scores (space-separated).

Example:
```
  A   C   D   E   F
A 2.0 0.3 -1.1 -2.7 -1.4
C 0.3 2.3  1.0 -1.1  0.2
...
```
The program will automatically convert this score matrix into a cost matrix (by subtracting each score from the maximum).
All amino acids in the matrix must be single uppercase letters and symmetric across rows and columns.

### Output Format

Results are saved in a `.tsv` file with at least the following columns:
```
query	match
```
If `v_gene` or `j_gene` information is present in the match entries, additional columns will be included:
```
query	match	v_gene	j_gene
```


## Installation

### Requirements
- **C++ Compiler:** C++17+
- **CMake:** Build system

### Build Instructions
```sh
git clone https://github.com/yourusername/TCRtrie.git
cd TCRtrie
mkdir build && cd build
cmake ..
make
```

## Contributing
If you encounter any bugs or have suggestions for improvements, please create an issue or submit a pull request on GitHub.
