# Constructors ---------------------------------------------------------------

#' Build Trie from AIRR TSV
#' @param path TSV with columns: junction_aa, v_call, j_call
#' @export
build_trie_airr <- function(path) {
  trie_create_from_airr(path)
}

#' Build Trie from vectors
#' @param seqs character vector of CDR3 AA
#' @param v_genes optional character vector (same length as seqs)
#' @param j_genes optional character vector (same length as seqs)
#' @export
build_trie_strings <- function(seqs, v_genes = NULL, j_genes = NULL) {
  seqs <- as.character(seqs)
  if (!is.null(v_genes)) v_genes <- as.character(v_genes)
  if (!is.null(j_genes)) j_genes <- as.character(j_genes)
  trie_create_from_strings(seqs, v_genes, j_genes)
}

# Config ---------------------------------------------------------------------

#' @export
trie_use_matrix <- function(trie, matrix_path, deletion_score = NULL) {
  trie_load_subst_matrix(trie, matrix_path)
  if (!is.null(deletion_score)) trie_set_deletion_score(trie, deletion_score)
  invisible(trie)
}

#' @export
trie_set_maxlen <- function(trie, n) {
  trie_set_max_query_length(trie, as.integer(n))
  invisible(trie)
}

#' @export
trie_print_matrix_df <- function(trie) {
  trie_print_matrix(trie)
  invisible(NULL)
}

# Searches -------------------------------------------------------------------

#' @export
trie_search_leven <- function(trie, query, max_edits) {
  trie_search(trie, as.character(query), as.integer(max_edits))
}

#' @export
trie_search_leven_any <- function(trie, query, max_edits) {
  trie_search_any(trie, as.character(query), as.integer(max_edits))
}

#' @export
trie_search_leven_batch <- function(trie, queries, max_edits) {
  trie_search_batch(trie, as.character(queries), as.integer(max_edits))
}

#' @export
trie_search_airr_df <- function(trie, query,
                                max_sub, max_ins, max_del,
                                max_edits = NULL,
                                v_filter = NULL, j_filter = NULL) {
  trie_search_airr(trie, as.character(query),
                   as.integer(max_sub), as.integer(max_ins), as.integer(max_del),
                   if (is.null(max_edits)) NULL else as.integer(max_edits),
                   v_filter, j_filter)
}

#' @export
trie_search_matrix_df <- function(trie, query, max_cost,
                                  v_filter = NULL, j_filter = NULL) {
  trie_search_matrix(trie, as.character(query), as.numeric(max_cost),
                     v_filter, j_filter)
}

#' @export
trie_search_for_all_airr_df <- function(trie, queries,
                                        max_sub, max_ins, max_del,
                                        max_edits = NULL,
                                        v_filter = NULL, j_filter = NULL) {
  trie_search_for_all_airr(trie, as.character(queries),
                           as.integer(max_sub), as.integer(max_ins), as.integer(max_del),
                           if (is.null(max_edits)) NULL else as.integer(max_edits),
                           v_filter, j_filter)
}

#' @export
trie_search_for_all_matrix_df <- function(trie, queries, max_cost,
                                          v_filter = NULL, j_filter = NULL) {
  trie_search_for_all_matrix(trie, as.character(queries), as.numeric(max_cost),
                             v_filter, j_filter)
}

# Cluster coverage -----------------------------------------------------------

#' @export
trie_cluster_usage_df <- function(trie, cluster,
                                  max_sub, max_ins, max_del,
                                  max_edits = NULL,
                                  v_filter = NULL, j_filter = NULL) {
  trie_cluster_usage(trie, as.character(cluster),
                     as.integer(max_sub), as.integer(max_ins), as.integer(max_del),
                     if (is.null(max_edits)) NULL else as.integer(max_edits),
                     v_filter, j_filter)
}

#' @export
trie_cluster_usage_matrix_df <- function(trie, cluster, max_cost,
                                         v_filter = NULL, j_filter = NULL) {
  trie_cluster_usage_matrix(trie, as.character(cluster), as.numeric(max_cost),
                            v_filter, j_filter)
}
