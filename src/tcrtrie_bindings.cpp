#include <Rcpp.h>
#include "Trie.h"

using namespace Rcpp;


static DataFrame to_df(const std::vector<AIRREntity>& v) {
    size_t n = v.size();
    CharacterVector cdr3(n), v_call(n), j_call(n);
    NumericVector dist(n);
    for (size_t i = 0; i < n; ++i) {
        cdr3[i] = v[i].junctionAA;
        v_call[i] = v[i].vGene;
        j_call[i] = v[i].jGene;
        dist[i] = v[i].distance;
    }
    return DataFrame::create(
            _["junction_aa"] = cdr3,
            _["v_call"]      = v_call,
            _["j_call"]      = j_call,
            _["distance"]    = dist,
            _["stringsAsFactors"] = false
    );
}

static DataFrame to_df_set(const std::unordered_set<AIRREntity>& s) {
    std::vector<AIRREntity> v;
    v.reserve(s.size());
    for (const auto& e : s) v.push_back(e);
    return to_df(v);
}

SEXP trie_create_from_airr(std::string data_path) {
    Rcpp::XPtr<Trie> ptr(new Trie(data_path), true);
    return ptr;
}

SEXP trie_create_from_strings(CharacterVector seqs,
                              Nullable<CharacterVector> v_genes = R_NilValue,
                              Nullable<CharacterVector> j_genes = R_NilValue) {
    std::vector<std::string> s; s.reserve(seqs.size());
    for (auto x : seqs) s.emplace_back(as<std::string>(x));

    std::vector<std::string> v, j;
    if (v_genes.isNotNull()) {
        CharacterVector vv(v_genes);
        if ((int)vv.size() != (int)s.size())
            stop("v_genes length must match seqs length");
        v.reserve(vv.size());
        for (auto x : vv) v.emplace_back(as<std::string>(x));
    } else {
        v.assign(s.size(), "");
    }

    if (j_genes.isNotNull()) {
        CharacterVector jj(j_genes);
        if ((int)jj.size() != (int)s.size())
            stop("j_genes length must match seqs length");
        j.reserve(jj.size());
        for (auto x : jj) j.emplace_back(as<std::string>(x));
    } else {
        j.assign(s.size(), "");
    }

    Rcpp::XPtr<Trie> ptr(new Trie(s, v, j), true);
    return ptr;
}

void trie_load_subst_matrix(SEXP xp, std::string matrix_path) {
    XPtr<Trie> ptr(xp);
    ptr->LoadSubstitutionMatrix(matrix_path);
}

void trie_set_deletion_score(SEXP xp, double score) {
    XPtr<Trie> ptr(xp);
    ptr->SetDeletionScore(static_cast<float>(score));
}

void trie_set_max_query_length(SEXP xp, int n) {
    XPtr<Trie> ptr(xp);
    ptr->SetMaxQueryLength(n);
}

void trie_print_matrix(SEXP xp) {
    XPtr<Trie> ptr(xp);
    ptr->PrintMatrix();
}

CharacterVector trie_search(SEXP xp, std::string query, int max_edits) {
    XPtr<Trie> ptr(xp);
    auto res = ptr->Search(query, max_edits);
    return wrap(res);
}

bool trie_search_any(SEXP xp, std::string query, int max_edits) {
    XPtr<Trie> ptr(xp);
    return ptr->SearchAny(query, max_edits);
}

List trie_search_batch(SEXP xp, CharacterVector queries, int max_edits) {
    XPtr<Trie> ptr(xp);
    std::vector<std::string> qs; qs.reserve(queries.size());
    for (auto q : queries) qs.emplace_back(as<std::string>(q));
    auto m = ptr->Search(qs, max_edits);

    List out(m.size());
    CharacterVector names(m.size());
    size_t i = 0;
    for (auto &kv : m) {
        names[i] = kv.first;
        out[i] = wrap(kv.second);
        ++i;
    }
    out.attr("names") = names;
    return out;
}

DataFrame trie_search_airr(SEXP xp, std::string query,
                           int max_sub, int max_ins, int max_del,
                           Nullable<int> max_edits = R_NilValue,
                           Nullable<std::string> v_filter = R_NilValue,
                           Nullable<std::string> j_filter = R_NilValue) {
    XPtr<Trie> ptr(xp);
    std::optional<int> me;
    if (max_edits.isNotNull()) me = as<int>(max_edits);
    std::optional<std::string> vf, jf;
    if (v_filter.isNotNull()) vf = as<std::string>(v_filter);
    if (j_filter.isNotNull()) jf = as<std::string>(j_filter);

    auto res = ptr->SearchAIRR(query, max_sub, max_ins, max_del, me, vf, jf);
    return to_df(res);
}

DataFrame trie_search_matrix(SEXP xp, std::string query, double max_cost,
                             Nullable<std::string> v_filter = R_NilValue,
                             Nullable<std::string> j_filter = R_NilValue) {
    XPtr<Trie> ptr(xp);
    std::optional<std::string> vf, jf;
    if (v_filter.isNotNull()) vf = as<std::string>(v_filter);
    if (j_filter.isNotNull()) jf = as<std::string>(j_filter);
    auto res = ptr->SearchWithMatrix(query, static_cast<float>(max_cost), vf, jf);
    return to_df(res);
}

List trie_search_for_all_airr(SEXP xp, CharacterVector queries,
                              int max_sub, int max_ins, int max_del,
                              Nullable<int> max_edits = R_NilValue,
                              Nullable<std::string> v_filter = R_NilValue,
                              Nullable<std::string> j_filter = R_NilValue) {
    XPtr<Trie> ptr(xp);
    std::vector<std::string> qs; qs.reserve(queries.size());
    for (auto q : queries) qs.emplace_back(as<std::string>(q));

    std::optional<int> me;
    if (max_edits.isNotNull()) me = as<int>(max_edits);
    std::optional<std::string> vf, jf;
    if (v_filter.isNotNull()) vf = as<std::string>(v_filter);
    if (j_filter.isNotNull()) jf = as<std::string>(j_filter);

    auto m = ptr->SearchForAll(qs, max_sub, max_ins, max_del, me, vf, jf);

    List out(m.size());
    CharacterVector names(m.size());
    size_t i = 0;
    for (auto &kv : m) {
        names[i] = kv.first;
        out[i] = to_df(kv.second);
        ++i;
    }
    out.attr("names") = names;
    return out;
}

List trie_search_for_all_matrix(SEXP xp, CharacterVector queries, double max_cost,
                                Nullable<std::string> v_filter = R_NilValue,
                                Nullable<std::string> j_filter = R_NilValue) {
    XPtr<Trie> ptr(xp);
    std::vector<std::string> qs; qs.reserve(queries.size());
    for (auto q : queries) qs.emplace_back(as<std::string>(q));

    std::optional<std::string> vf, jf;
    if (v_filter.isNotNull()) vf = as<std::string>(v_filter);
    if (j_filter.isNotNull()) jf = as<std::string>(j_filter);

    auto m = ptr->SearchForAllWithMatrix(qs, static_cast<float>(max_cost), vf, jf);

    List out(m.size());
    CharacterVector names(m.size());
    size_t i = 0;
    for (auto &kv : m) {
        names[i] = kv.first;
        out[i] = to_df(kv.second);
        ++i;
    }
    out.attr("names") = names;
    return out;
}

DataFrame trie_cluster_usage(SEXP xp, CharacterVector cluster,
                             int max_sub, int max_ins, int max_del,
                             Nullable<int> max_edits = R_NilValue,
                             Nullable<std::string> v_filter = R_NilValue,
                             Nullable<std::string> j_filter = R_NilValue) {
    XPtr<Trie> ptr(xp);
    std::vector<std::string> cl; cl.reserve(cluster.size());
    for (auto q : cluster) cl.emplace_back(as<std::string>(q));

    std::optional<int> me;
    if (max_edits.isNotNull()) me = as<int>(max_edits);
    std::optional<std::string> vf, jf;
    if (v_filter.isNotNull()) vf = as<std::string>(v_filter);
    if (j_filter.isNotNull()) jf = as<std::string>(j_filter);

    auto res = ptr->ClusterUsage(cl, max_sub, max_ins, max_del, me, vf, jf);
    return to_df_set(res);
}

DataFrame trie_cluster_usage_matrix(SEXP xp, CharacterVector cluster,
                                    double max_cost,
                                    Nullable<std::string> v_filter = R_NilValue,
                                    Nullable<std::string> j_filter = R_NilValue) {
    XPtr<Trie> ptr(xp);
    std::vector<std::string> cl; cl.reserve(cluster.size());
    for (auto q : cluster) cl.emplace_back(as<std::string>(q));

    std::optional<std::string> vf, jf;
    if (v_filter.isNotNull()) vf = as<std::string>(v_filter);
    if (j_filter.isNotNull()) jf = as<std::string>(j_filter);

    auto res = ptr->ClusterUsageWithMatrix(cl, static_cast<float>(max_cost), vf, jf);
    return to_df_set(res);
}
