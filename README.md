# TCRtrie

TCRtrie is a tool for approximate search in TCR repertoires based on a CDR3 index.
It can be used both for searching user-provided repertoires and for searching the
[VDJdb](https://vdjdb.com) database.

The library supports two search modes:

- edit-distance-based search with bounded substitutions, insertions, and deletions;
- matrix-based search where substitutions are scored using an amino acid substitution matrix.

The core idea is to build a trie index over CDR3 amino acid sequences and use it for fast approximate matching.

## Main objects

### `VDJdb`

`VDJdb` is a **lazy object**.
It is initialized on first use and is **not callable**.
Use it like this:

```python
from tcrtrie import VDJdb

VDJdb.search(...)
```

Do **not** write:

```python
VDJdb()
```

### `Trie`

`Trie` is the low-level class for building an index from your own AIRR-like TSV file
or from in-memory sequence arrays.

## Examples

### Search in VDJdb with edit distance

```python
from tcrtrie import VDJdb

df = VDJdb.search(  # returns pandas.DataFrame
    query="CASSEGTDGYTF",
    maxSubstitution=2,
    maxInsertion=1,
    maxDeletion=1,
    maxEdits=2,
    vGeneFilter="TRBV19*01",
    jGeneFilter="TRBJ1-2*01",
    numThreads=8,
    detailed=True,
)

df
```

### Batch search in VDJdb

```python
from tcrtrie import VDJdb

VDJdb.searchForAll(
    queries=["CASSEGTDGYTF", "CAISTGDSNQPQHF"],
    maxSubstitution=2,
    maxInsertion=1,
    maxDeletion=1,
    maxEdits=2,
    vGeneFilters=["TRBV19*01", "TRBV6-6*01"],
    jGeneFilters=["TRBJ1-2*01", "TRBJ1-5*01"],
    numThreads=8,
    detailed=True,
)
```

### Search in VDJdb with a substitution matrix

```python
from tcrtrie import VDJdb

VDJdb.searchWithMatrix(
    query="CASSEGTDGYTF",
    maxCost=12,
    detailed=True,
)
```

### Build a trie from your own TSV file

```python
from tcrtrie import Trie

trie = Trie("my_repertoire.tsv")

hits = trie.SearchIndices(
    query="CASSEGTDGYTF",
    maxSubstitution=2,
    maxInsertion=1,
    maxDeletion=1,
    maxEdits=2,
)

hits
```

### Load a custom substitution matrix

```python
from tcrtrie import Trie

trie = Trie("my_repertoire.tsv")
trie.LoadSubstitutionMatrix("my_matrix.txt", delimiter="", gapFactor=1.5)

hits = trie.SearchIndicesWithMatrix(
    query="CASSEGTDGYTF",
    maxCost=12,
)

hits
```

## Installation

### Requirements

- Python 3.8+
- `pip`
- C++17-compatible compiler
- CMake

On Windows, you may need to install Microsoft C++ [Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/) if no suitable compiler is available.
### Install from PyPI

```bash
pip install tcrtrie
```

### Install directly from GitHub

```bash
python -m pip install git+https://github.com/MikePodsytnik/TCRtrie@0.2.0-tcrtriepy
```

### Clone and build locally

```bash
git clone --branch TCRtriePy https://github.com/MikePodsytnik/TCRtrie.git
cd TCRtrie
python -m pip install --upgrade pip setuptools wheel scikit-build-core pybind11
python -m pip install .
```

## VDJdb database management

TCRtriePy does not update VDJdb automatically in order to avoid silent changes in scientific results.
Database updates are performed explicitly via a command-line tool.

The `tcrtrie-vdjdb-update` command downloads and installs a selected VDJdb release into the local cache
(`~/.cache/tcrtrie/vdjdb`). The cached version is then used by the `VDJdb` object in Python.

Install the latest available VDJdb release:
```{bash}
tcrtrie-vdjdb-update
```

Install a specific VDJdb version (recommended for reproducibility):
```{bash}
tcrtrie-vdjdb-update --tag 2025-12-29
```

List available VDJdb releases:
```{bash}
tcrtrie-vdjdb-update --list | head -n 10
```

After updating the database, restart the Python process or Jupyter kernel
to ensure the new version is used.

## Input data format

### Repertoire TSV format

`Trie` expects an AIRR-like tab-separated file.

Required column:

- `junction_aa`

Optional columns:

- `v_call`
- `j_call`
- `__group_id`

Minimal example:

```tsv
junction_aa	v_call	j_call
CASSEGTDGYTF	TRBV19*01	TRBJ1-2*01
CAISTGDSNQPQHF	TRBV6-6*01	TRBJ1-5*01
```

Notes:

- the file must be tab-separated;
- amino acid sequences are read from `junction_aa`;
- Sequences must be uppercase and contain only the 20 standard amino acid symbols (`ACDEFGHIKLMNPQRSTVWY`).
- `v_call` and `j_call` are optional, but they are required if you want to use V/J filtering;

### Substitution matrix format

The substitution matrix must be a square matrix over amino acid symbols.

Requirements:

- all 20 standard amino acids must be present;
- the matrix may contain either 20 labels or 21 labels if gap (`-`) is provided explicitly;
- row and column labels must match;
- diagonal values must be strictly greater than every other value in the corresponding row and column.

Whitespace-separated matrices are supported by default.
If your matrix uses another separator, pass it via the `delimiter` argument.

Example:

```text
   A  R  N  D ...
A  4 -1 -2 -2 
R -1  5  0 -2 
N -2  0  6  1 
D -2 -2  1  6 
...
```

#### How gap scores are handled

If the matrix already contains `-`, those values are used directly.

If a gap column is not provided, the gap score for amino acid `a` is derived as the negative diagonal score, i.e. `gap(a) = -score(a, a) * gapFactor`. 
Then the same synthesized value is written to both `aa -> -` and `- -> aa`.
By default, `gapFactor=1.0`.

After that, the score matrix is converted into an internal non-negative cost matrix used by matrix-based search.
