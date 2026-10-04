# TCRtrie

`TCRtrie` is a command-line tool for fast approximate matching of T-cell receptor (TCR) repertoires using trie-based sequence search.

It compares two clonotype tables in TSV format:
- the first table is treated as the query repertoire;
- the second table is treated as the target repertoire;
- for each query clonotype, all target clonotypes within the specified search radius are reported.

The result file contains information about both the query clonotype and the matched target clonotype. The same target clonotype may therefore appear multiple times if it matches multiple queries.

## Quick start

Search for sequences differing by at most one substitution:

```bash
tcrtrie query.tsv target.tsv --max-sub 1 --max-edits 1
```

Write results to a custom file:

```bash
tcrtrie query.tsv target.tsv \
  --out matches.tsv \
  --max-sub 1 \
  --match-v \
  --match-j
```

Show all available options:
```
tcrtrie --help
```

## CLI arguments

| Flag             |        Type | Default                | Description                                                                                                |
|------------------|------------:|------------------------|------------------------------------------------------------------------------------------------------------|
| `query_tsv`      |  positional | —                      | First input clonotype table (query repertoire, TSV).                                                       |
| `target_tsv`     |  positional | —                      | Second input clonotype table (target repertoire, TSV).                                                     |
| `--out`          |       `str` | `match_result.tsv`     | Output TSV path.                                                                                           |
| `--max-sub`      |       `int` | `1`                    | Maximum substitutions allowed.                                                                             |
| `--max-ins`      |       `int` | `0`                    | Maximum insertions allowed.                                                                                |
| `--max-del`      |       `int` | `0`                    | Maximum deletions allowed.                                                                                 |
| `--max-edits`    |       `int` | `max-sub + max-ins + max-del` | Maximum total number of edit operations allowed.                                                           |
| `--matrix-path`  |       `str` | —                      | Path to a substitution matrix file. If set, matrix search mode is used.                                    |
| `--matrix-v` | `str` | — | Path to substitution matrix for the V-derived CDR3 region. |
| `--matrix-ndn` | `str` | — | Path to substitution matrix for the NDN-derived CDR3 region. |
| `--matrix-j` | `str` | — | Path to substitution matrix for the J-derived CDR3 region. |
| `--max-cost`     | `int/float` | `6`                    | Maximum total alignment cost allowed in matrix mode.                                                       |
| `--match-v`      |      `flag` | off                    | Require V gene to match when counting a match.                                                             |
| `--match-j`      |      `flag` | off                    | Require J gene to match when counting a match.                                                             |
| `--gene`         |       `str` | —                      | Chain filter applied while reading input tables.                                                           |
| `--species`      |       `str` | —                      | Species filter applied while reading input tables.                                                         |
| `--epitope`      |       `str` | —                      | Epitope filter applied while reading input tables. |
| `--threads`      |       `int` | `4`                    | Maximum number of worker threads for search.                                                               |
| `--junction-col` |       `str` | `junction_aa or cdr3`  | Column name for the CDR3 amino-acid sequence. If not specified, the reader tries common alternatives.      |
| `--v-end-col` |       `str` | `v.end or vEnd`  |       |
| `--j-start-col` |       `str` | `j.start or jStart`  |       |
| `--cdr3fix-col` |       `str` | `cdr3fix`  |       |
| `--v-col`        |       `str` | `v_call or v.segm`     | Column name for the V gene call. If not specified, the reader tries common alternatives.                   |
| `--j-col`        |       `str` | `j_call or j.segm`     | Column name for the J gene call. If not specified, the reader tries common alternatives.                   |
| `--epitope-col`  |       `str` | `antigen.epitope`      | Column name for the epitope. This column is optional.                                                      |
| `--species-col`  |       `str` | `species`              | Column name for the species label.                                                                         |
| `--chain-col`    |       `str` | `gene`                 | Column name for the chain label, for example `TRB` or `TRA`.                                               |
| `--align`        |      `flag` | off                    | Add alignment information to the output.                                                                   |

## Search modes

`tcrtrie` supports three search modes.

### 1. Bounded edit mode

This is the default mode if `--matrix-path` is not provided.

The search is controlled by:
- `--max-sub`
- `--max-ins`
- `--max-del`
- `--max-edits`

Example:

```bash
tcrtrie query.tsv target.tsv \
  --out match_result.tsv \
  --max-sub 1 \
  --max-ins 0 \
  --max-del 0 \
  --max-edits 1 \
  --threads 4 \
  --match-v \
  --match-j
```

### 2. Substitution-matrix mode

If `--matrix-path` is provided, the tool switches to matrix-based matching.

The search is controlled by:
- `--matrix-path`
- `--max-cost`

Example:

```bash
tcrtrie query.tsv target.tsv \
  --matrix-path blosum62.txt \
  --max-cost 6 \
  --threads 4 \
  --align
```

### 3. Region-specific substitution-matrix mode

Region-aware search uses separate substitution matrices for the
V-derived, N(D)N-derived, and J-derived regions of CDR3.

All three matrices must be provided:

```bash
tcrtrie query.tsv target.tsv \
  --matrix-v matrix_v.txt \
  --matrix-ndn matrix_ndn.txt \
  --matrix-j matrix_j.txt \
  --max-cost 6
```

Regional boundaries are read from `v.end` and `j.start` columns by default. They can also be supplied through the `cdr3fix` annotation or configured using `--v-end-col`, `--j-start-col`, and `--cdr3fix-col`.

### Search against VDJdb
If you want to search against VDJdb, you can use it directly:
- Download a release from: https://github.com/antigenomics/vdjdb-db/releases
- Unpack the archive
- Use vdjdb.txt as your target file (no need to modify columns)

Example of a 1-sub search against VDJdb:
```bash
tcrtrie query.tsv vdjdb.txt 
```

## Installation

#### Requirements
- Python 3.8+
- pip
- CMake >= 3.16
- C++17-compatible compiler

`tcrtrie-cli` is distributed as source code and compiled locally during installation.
A working C++17 compiler must therefore be available on your system.

### 1. Install from PyPI

```bash
pip install tcrtrie-cli
```

### 2. Install directly from GitHub


```bash
pip install "git+https://github.com/MikePodsytnik/TCRtrie.git"
```

### 3. Clone the repository and install locally

```bash
git clone https://github.com/MikePodsytnik/TCRtrie.git
cd TCRtrie
pip install .
```

After installation, the CLI should be available as:

```bash
tcrtrie --help
```

If installation fails, check that:
- `cmake` is available in your shell
- a working C++ compiler is installed
- Python and `pip` are available in the environment where you run the installation


## Input and output format

`tcrtrie` takes **two input TSV files**: a **query repertoire** and a **target repertoire**.  
Each row is a single clonotype. The tool iterates over clonotypes from the first file and searches for matches in the second one.

Both files must be **tab-separated** and contain a **header row**.

### Input columns

Only the CDR3 amino-acid sequence column is required. By default, `tcrtrie` tries to detect it as:

- `junction_aa`
- `cdr3`

You can also set it explicitly with `--junction-col`.

Other columns are optional and are used only when needed:

- V gene: `--v-col` (*v_call* or *v.segm*)
- J gene: `--j-col` (*j_call* or *j.segm*)
- epitope: `--epitope-col` (*antigen.epitope*)
- species: `--species-col` (*species*)
- chain: `--chain-col` (*gene*)

If a filter is requested, the corresponding column must be present:

- `--match-v` → V column
- `--match-j` → J column
- `--epitope` → epitope column
- `--species` → species column
- `--gene` → chain column

Dataset-level filters (`--gene`, `--species`, `--epitope`) are applied **before** trie construction and **before** matching.

### Output file

The output is a TSV where **each row is one query–target match**.  
If several target clonotypes fall into the allowed radius for one query clonotype, all of them are written as separate rows. The same target clonotype may therefore appear multiple times.

The output contains fields from both the query and target clonotypes, plus match metadata.

Depending on the mode, it also includes:

- **edit mode**: distance, substitutions, insertions, deletions
- **matrix mode**: match cost
- **optional alignment mode**: query alignment, target alignment

So the output is **not** a deduplicated repertoire intersection, but a full list of all matched query–target pairs.

## Citation

If you use TCRtrie in your research, please cite:

> M. Podsytnik, E. Vlasova, and M. Shugay, "TCRtrie: Scalable Fuzzy Retrieval with Configurable Similarity for T-Cell Receptor Repertoires," 2026.

## License

TCRtrie is released under the [MIT License](LICENSE).

