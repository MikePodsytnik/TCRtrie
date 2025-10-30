# TCRtrieR: Fast R Bindings for TCR Sequence Matching with Trie

`tcrtrieR` provides R bindings for the C++ **TCRtrie** library, enabling fast approximate matching and similarity search over T‑cell receptor (TCR) CDR3 sequences using a trie (prefix tree) data structure. Searches can be performed with Levenshtein-style edit limits or with substitution matrices (e.g., BLOSUM), and can be filtered by V/J genes.

---

## Overview

`tcrtrieR` helps analyze TCR repertoires by:
- Building a trie index from an AIRR TSV or from vectors of CDR3 sequences (optionally with V and J metadata).
- Performing approximate matching with per‑operation limits (substitutions/insertions/deletions).
- Using substitution matrices with configurable deletion (gap) penalties.
- Filtering by **V** and **J** gene segments.
- Running batch searches with parallel execution underneath.

---

## Installation

### Requirements
- **R ≥ 4.2**
- A C++17 toolchain (clang or gcc)
- R packages: `Rcpp`, `devtools` (installed automatically if missing)

### Install from GitHub

```r
install.packages("remotes")
remotes::install_github("MikePodsytnik/TCRtrie", ref = "0.1.0-tcrtrieR")

```

### Install from cloned sources

```bash
git clone https://github.com/MikePodsytnik/tcrtrieR.git
cd tcrtrieR

R -q -e "Rcpp::compileAttributes('.')"

R CMD INSTALL .
```

---

## R API Usage

```r
library(tcrtrieR)

# Build Trie from AIRR TSV (expects columns: junction_aa, v_call, j_call)
tr <- build_trie_airr("vdjdb_airr.tsv")

# Levenshtein-based AIRR search with per-operation limits and V/J filters
airr_results <- trie_search_airr_df(
  tr,
  query    = "CASSEGTDGYTF",
  max_sub  = 2,
  max_ins  = 1,
  max_del  = 1,
  v_filter = "TRBV4-1*01",
  j_filter = "TRBJ1-2*01"
)

# Print results (junction AA, V, J, distance)
for (i in seq_len(nrow(airr_results))) {
  cat(
    airr_results$junction_aa[i], " ",
    airr_results$v_call[i], " ",
    airr_results$j_call[i], " ",
    airr_results$distance[i], "
"
  )
}
```

### Example: Substitution Matrix

```r
library(tcrtrieR)

tr <- build_trie_airr("vdjdb_airr.tsv")

# Load a BLOSUM-like matrix and configure deletion score if the matrix lacks a gap column
trie_use_matrix(tr, "blosum.txt", deletion_score = -5)

matrix_results <- trie_search_matrix_df(
  tr,
  query    = "CASSLATDGYTF",
  max_cost = 5.0,
  v_filter = "TRBV5-6*01",
  j_filter = "TRBJ1-2*01"
)

head(matrix_results)
```

---

## Key Functions

- `build_trie_airr(path)`  
  Build a Trie index from an AIRR TSV file (columns: `junction_aa`, `v_call`, `j_call`).

- `build_trie_strings(seqs, v_genes = NULL, j_genes = NULL)`  
  Build from vectors; `v_genes` and `j_genes` are optional vectors of the same length as `seqs`.

- `trie_search_leven(trie, query, max_edits)`  
  Simple Levenshtein search returning matching CDR3 strings.

- `trie_search_airr_df(trie, query, max_sub, max_ins, max_del, max_edits = NULL, v_filter = NULL, j_filter = NULL)`  
  AIRR-aware search returning a `data.frame` with `junction_aa`, `v_call`, `j_call`, `distance`.

- `trie_search_matrix_df(trie, query, max_cost, v_filter = NULL, j_filter = NULL)`  
  Search using a substitution matrix (cost threshold). Returns a `data.frame`.

- `trie_search_for_all_airr_df(trie, queries, ...)`  
  Batch AIRR search for multiple queries; returns a named list of `data.frame`s.

- `trie_search_for_all_matrix_df(trie, queries, ...)`  
  Batch search with a substitution matrix; returns a named list of `data.frame`s.

- `trie_cluster_usage_df(trie, cluster, max_sub, max_ins, max_del, ...)`  
  Return all AIRR records matched to the provided cluster (as a `data.frame`).

- `trie_use_matrix(trie, matrix_path, deletion_score = NULL)`  
  Load a substitution matrix and optionally set a deletion (gap) penalty used in matrix-based search.

- `trie_set_maxlen(trie, n)`  
  Set the maximum allowed query length (default is 32).

---

## Substitution Matrix Format

Matrix files should follow a BLOSUM-like layout, optionally including a first row/column for gaps `-`:

```
    -   A    C    D    E
-  -3 -0.8 -0.7 -1.3 -0.1
A   2.0 0.3 -1.1 -2.7 -1.4
C   0.3 2.3  1.0 -1.1  0.2
...
```

If the matrix does not define a gap row/column, use `trie_use_matrix(tr, path, deletion_score = <value>)` to configure the deletion (gap) cost.

---

## Notes

- Default `maxQueryLength` is **32**; increase it via `trie_set_maxlen(tr, new_limit)` if needed.
- Batch functions (`trie_search_for_all_*`) utilize parallel execution via `std::async` internally.
- The package is implemented in C++17 and bound to R via [Rcpp](https://cran.r-project.org/package=Rcpp).

---

## Troubleshooting

- `Rcpp.h: No such file or directory` → Install `Rcpp` (`install.packages("Rcpp")`) and reinstall the package.
- `Illegal instruction` on execution → Reinstall without CPU‑specific flags (ensure `TCRTRIER_FASTMATH` is unset).
- Windows toolchain issues → Install **Rtools** for your R version and verify with `pkgbuild::has_build_tools()`.

---

## Contributing

Issues and pull requests are welcome at: https://github.com/MikePodsytnik/tcrtrieR

---

