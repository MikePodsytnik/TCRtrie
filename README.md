# TCRtriePy: Efficient Python Bindings for TCR Sequence Matching with Trie

TCRtriePy provides fast Python bindings for the TCRtrie C++ library, leveraging a trie (prefix tree) data structure for approximate matching and efficient similarity searches within T-cell receptor (TCR) sequence datasets and [pybind11](https://github.com/pybind/pybind11) for seamless Python integration.
## Overview

TCRtriePy addresses the complexity of comparing TCR sequences, especially CDR3 regions, by:

* Constructing a trie from a collection of CDR3 sequences.
* Performing approximate matches using Levenshtein-based distances or substitution matrices.
* Accelerating searches through multithreading.
* Supporting optional filtering by V and J gene segments.

## Installation

### Requirements

* Python ≥ 3.8
* C++ Compiler supporting C++17
* CMake ≥ 3.16

### Install from Source

```bash
git clone --branch TCRtriePy https://github.com/MikePodsytnik/TCRtrie.git
cd TCRtrie
pip install --upgrade pip setuptools wheel scikit-build-core pybind11
pip install .
```

## Python API Usage

```python
from tcrtrie import Trie

# Initialize Trie from AIRR file
trie = Trie("vdjdb_airr.tsv")

# Advanced AIRR-based search
airr_results = trie.SearchAIRR(
   query="CASSLGTDGYTF",
   maxSubstitution=2,
   maxInsertion=1,
   maxDeletion=1,
   vGeneFilter="TRBV7-9",
   jGeneFilter="TRBJ2-1"
)

for r in airr_results:
   print(r.junctionAA, r.vGene, r.jGene, r.distance)
```

## Key Features

`Search(query, maxEdits)`

Finds sequences matching a query within a specified edit distance.

`SearchAIRR(query, maxSubstitution, maxInsertion, maxDeletion, vGeneFilter, jGeneFilter)`

Advanced AIRR-compliant search with individual edit type constraints and optional gene filters.

`SearchWithMatrix(query, maxCost, vGeneFilter, jGeneFilter)`

Search using a substitution matrix with a defined cost threshold.

`SearchAny(query, maxEdits)`

Checks if any sequence matches the query within a defined edit distance.

`SearchForAll(queries, maxSubstitution, maxInsertion, maxDeletion, vGeneFilter, jGeneFilter)`

Performs batch searches with multithreading.

`LoadSubstitutionMatrix(matrixPath)`

Loads and converts a substitution matrix for use in searches.

`SetMaxQueryLength(newMaxQueryLength)`

Adjusts maximum allowed query length (default is 32).

## Substitution Matrix Format

Matrix files should follow a BLOSUM-like format:

```
  A   C   D   E   F
A 2.0 0.3 -1.1 -2.7 -1.4
C 0.3 2.3  1.0 -1.1  0.2
...
```

## Contributing

Feel free to submit issues or pull requests via GitHub if you encounter any problems or have suggestions for improvement.
