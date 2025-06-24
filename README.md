# TCRtriePy: Efficient Python Bindings for TCR Sequence Matching with Trie

TCRtriePy provides fast Python bindings for the TCRtrie C++ library, leveraging a trie (prefix tree) data structure for approximate matching and efficient similarity searches within T-cell receptor (TCR) sequence datasets and [pybind11](https://github.com/pybind/pybind11) for seamless Python integration.
## Overview

TCRtriePy addresses the complexity of comparing TCR sequences, especially CDR3 regions, by:

* Constructing a trie from a collection of CDR3 sequences.
* Performing approximate matches using Levenshtein-based distances or substitution matrices.
* Accelerating searches through multithreading.
* Supporting optional filtering by V and J gene segments.

## Installation

### Installation procedure
One can simply install the software out-of-the-box using [pip](https://pypi.org/project/pip/)

```bash
conda create -n tcrtrie
conda activate tcrtrie
pip install git+https://github.com/MikePodsytnik/TCRtrie@0.1.0-tcrtriepy
```

Or, in case of package version problems or other issues, clone the repository manually via git, create
corresponding [conda](https://docs.conda.io/en/latest/) environment and install directly from sources:

```{bash}
git clone --branch TCRtriePy https://github.com/MikePodsytnik/TCRtrie.git
cd TCRtrie
conda create -n tcrtrie
conda activate tcrtrie
pip install --upgrade pip setuptools wheel scikit-build-core pybind11
pip install .
```

> For this method, ensure your system has a C++ Compiler supporting C++17 and CMake ≥ 3.16 installed


## Python API Usage

```python
from tcrtrie import Trie

# Initialize Trie from AIRR file
trie = Trie("vdjdb_airr.tsv")

# Levenshtein distance AIRR-based search
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
Example using Substitution Matrix and [VDJdb](https://vdjdb.cdr3.net/)
```python
from tcrtrie import VDJdb

# If there is no corresponding column in the matrix.txt
VDJdb.SetDeletionScore(-5)
VDJdb.LoadSubstitutionMatrix("../blosum.txt")

matrix_results = VDJdb.SearchWithMatrix(
  query="CASSLATDGYTF",
  maxCost=5.0,
  vGeneFilter="TRBV5-6*01",
  jGeneFilter="TRBJ1-2*01"
)

for r in matrix_results:
  print(r.junctionAA, r.vGene, r.jGene, "%.3f" %(r.distance))
```

## Key Features

* `Search(query, maxEdits)`

    Finds sequences matching a query within a specified edit distance.


* `SearchAIRR(query, maxSubstitution, maxInsertion, maxDeletion, vGeneFilter, jGeneFilter)`

    Advanced AIRR-compliant search with individual edit type constraints and optional gene filters.


* `SearchWithMatrix(query, maxCost, vGeneFilter, jGeneFilter)`

    Search using a substitution matrix with a defined cost threshold.


* `SearchAny(query, maxEdits)`

    Checks if any sequence matches the query within a defined edit distance.


* `SearchForAll(queries, maxSubstitution, maxInsertion, maxDeletion, vGeneFilter, jGeneFilter)`

    Performs batch searches with multithreading.


* `LoadSubstitutionMatrix(matrixPath)`

    Loads a substitution matrix and converts into substitution-cost matrix for use in searches.


* `SetDeletionScore(float deletionScore);`

    Adjusts deletion score an amino acids for use in SearchWithMatrix if you don't have it in matrix (default is -6).


* `SetMaxQueryLength(newMaxQueryLength)`

    Adjusts maximum allowed query length (default is 32).

## Substitution Matrix Format

Matrix files should follow a BLOSUM-like format, but may also have a column for the gap score.

```
    -   A    C    D    E
-  -3 -0.8 -0.7 -1.3 -0.1
A  2.0 0.3 -1.1 -2.7 -1.4
C  0.3 2.3  1.0 -1.1  0.2
...
```

## Contributing

Feel free to submit issues or pull requests via GitHub if you encounter any problems or have suggestions for improvement.
