#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "AirrParser.h"
#include "Trie.h"

namespace py = pybind11;

template <typename F, typename Obj, typename... Args>
auto call_without_gil(F&& f, Obj&& obj, Args&&... args) {
    py::gil_scoped_release release;
    return (std::forward<Obj>(obj)->*f)(std::forward<Args>(args)...);
}

PYBIND11_MODULE(_tcrtrie, m) {
    m.doc() = "Python bindings for TCRtrie C++ library";

    py::class_<AIRREntity>(
        m,
        "AIRREntity",
        R"pbdoc(
AIRR-style search result.

Attributes
----------
junctionAA : str
    Matched amino acid CDR3 sequence.
vGene : str
    V gene associated with the matched sequence.
jGene : str
    J gene associated with the matched sequence.
distance : int
    Edit distance between the query and the matched sequence.
)pbdoc")
        .def_readonly("junctionAA", &AIRREntity::junctionAA, "Matched amino acid CDR3 sequence.")
        .def_readonly("vGene", &AIRREntity::vGene, "V gene associated with the matched sequence.")
        .def_readonly("jGene", &AIRREntity::jGene, "J gene associated with the matched sequence.")
        .def_readonly("distance", &AIRREntity::distance, "Edit distance between query and target.")
        .def("__repr__", [](const AIRREntity& e) {
            return "<AIRREntity junctionAA='" + e.junctionAA + "' vGene='" + e.vGene +
                   "' jGene='" + e.jGene + "' distance=" + std::to_string(e.distance) + ">";
        });

    py::enum_<Trie::AlignmentOpType>(
        m,
        "AlignmentOpType",
        R"pbdoc(
Type of operation in an alignment.
)pbdoc")
        .value("Match", Trie::AlignmentOpType::Match)
        .value("Substitution", Trie::AlignmentOpType::Substitution)
        .value("Insertion", Trie::AlignmentOpType::Insertion)
        .value("Deletion", Trie::AlignmentOpType::Deletion)
        .export_values();

    py::class_<Trie::AlignmentOp>(
        m,
        "AlignmentOp",
        R"pbdoc(
Single alignment operation.

Attributes
----------
type : AlignmentOpType
    Operation type.
queryPos : int
    Position in the original query sequence to which the operation относится.
queryChar : str
    Character from the query sequence.
targetChar : str
    Character from the target sequence.
)pbdoc")
        .def_readonly("type", &Trie::AlignmentOp::type, "Operation type.")
        .def_readonly("queryPos", &Trie::AlignmentOp::queryPos, "Position in the original query sequence.")
        .def_readonly("queryChar", &Trie::AlignmentOp::queryChar, "Character from the query sequence.")
        .def_readonly("targetChar", &Trie::AlignmentOp::targetChar, "Character from the target sequence.")
        .def("__repr__", [](const Trie::AlignmentOp& op) {
            return "<AlignmentOp type=" + py::str(py::cast(op.type)).cast<std::string>() +
                   " queryPos=" + std::to_string(op.queryPos) +
                   " queryChar='" + std::string(1, op.queryChar) +
                   "' targetChar='" + std::string(1, op.targetChar) + "'>";
        });

    py::class_<Trie::AlignmentResult>(
        m,
        "AlignmentResult",
        R"pbdoc(
Detailed alignment result.

Attributes
----------
queryAligned : str
    Query string after alignment, including gap characters if needed.
targetAligned : str
    Target string after alignment, including gap characters if needed.
substitutions : int
    Number of substitutions in the alignment.
insertions : int
    Number of insertions in the alignment.
deletions : int
    Number of deletions in the alignment.
distance : float
    Total alignment score or distance.
ops : list[AlignmentOp]
    Detailed sequence of alignment operations.
)pbdoc")
        .def_readonly("queryAligned", &Trie::AlignmentResult::queryAligned, "Aligned query string.")
        .def_readonly("targetAligned", &Trie::AlignmentResult::targetAligned, "Aligned target string.")
        .def_readonly("substitutions", &Trie::AlignmentResult::substitutions, "Number of substitutions.")
        .def_readonly("insertions", &Trie::AlignmentResult::insertions, "Number of insertions.")
        .def_readonly("deletions", &Trie::AlignmentResult::deletions, "Number of deletions.")
        .def_readonly("distance", &Trie::AlignmentResult::distance, "Total alignment score or distance.")
        .def_readonly("ops", &Trie::AlignmentResult::ops, "Detailed alignment operations.")
        .def("__repr__", [](const Trie::AlignmentResult& r) {
            return "<AlignmentResult distance=" + std::to_string(r.distance) +
                   " substitutions=" + std::to_string(r.substitutions) +
                   " insertions=" + std::to_string(r.insertions) +
                   " deletions=" + std::to_string(r.deletions) +
                   " queryAligned='" + r.queryAligned +
                   "' targetAligned='" + r.targetAligned + "'>";
        });

    py::class_<Trie>(
        m,
        "Trie",
        R"pbdoc(
Trie-based index for TCR-like amino acid sequences.

The class supports:
- exact / bounded-edit search,
- batch search,
- search with optional V/J gene filters,
- substitution-matrix-based search,
- alignment reconstruction for query-target pairs or indexed hits.
)pbdoc")
        .def(py::init<>(), R"pbdoc(
Create an empty trie.
)pbdoc")
        .def(py::init<const std::string&>(),
             py::arg("dataPath"),
             R"pbdoc(
Create a trie from a data file.

Parameters
----------
dataPath : str
    Path to the input data file used to initialize the trie.
)pbdoc")
        .def(py::init<const Trie&>(),
             py::arg("other"),
             R"pbdoc(
Create a copy of another trie.

Parameters
----------
other : Trie
    Existing trie instance to copy.
)pbdoc")
        .def(py::init<const std::vector<std::string>&,
                      const std::vector<std::string>&,
                      const std::vector<std::string>&>(),
             py::arg("sequences"),
             py::arg("vGenes"),
             py::arg("jGenes"),
             R"pbdoc(
Create a trie from in-memory sequence and gene arrays.

Parameters
----------
sequences : list[str]
    Amino acid sequences to insert.
vGenes : list[str]
    V genes corresponding to sequences.
jGenes : list[str]
    J genes corresponding to sequences.
)pbdoc")

        .def("LoadSubstitutionMatrix",
             &Trie::LoadSubstitutionMatrix,
             py::arg("matrixPath"),
             py::arg("delimiter") = "",
             py::arg("gapFactor") = 1.5f,
             R"pbdoc(
Load a substitution matrix used by matrix-based search and alignment.

Parameters
----------
matrixPath : str
    Path to the matrix file.
delimiter : str, default ""
    Column delimiter. Empty string means whitespace-separated input.
gapFactor : float, default 1.5
    Factor used when synthesizing gap costs if they are not explicitly provided.
)pbdoc")

        .def("PrintMatrix",
             &Trie::PrintMatrix,
             R"pbdoc(
Print the currently loaded substitution matrix to standard output.
)pbdoc")

        .def("Search",
             [](Trie& self, const std::string& query, int maxEdits) {
                 using Method = std::vector<std::string> (Trie::*)(const std::string&, int);
                 return call_without_gil(static_cast<Method>(&Trie::Search), &self, query, maxEdits);
             },
             py::arg("query"),
             py::arg("maxEdits"),
             R"pbdoc(
Search sequences within a maximum edit distance.

Parameters
----------
query : str
    Query amino acid sequence.
maxEdits : int
    Maximum allowed edit distance.

Returns
-------
list[str]
    Matching sequences.
)pbdoc")

        .def("Search",
             [](Trie& self,
                const std::vector<std::string>& queries,
                int maxEdits,
                std::optional<std::size_t> numThreads) {
                 using Method = std::vector<std::vector<std::string>> (Trie::*)(
                     const std::vector<std::string>&,
                     int,
                     std::optional<std::size_t>);
                 return call_without_gil(static_cast<Method>(&Trie::Search),
                                         &self,
                                         queries,
                                         maxEdits,
                                         numThreads);
             },
             py::arg("queries"),
             py::arg("maxEdits"),
             py::arg("numThreads") = 4,
             R"pbdoc(
Search multiple queries within a maximum edit distance.

Parameters
----------
queries : list[str]
    Query amino acid sequences.
maxEdits : int
    Maximum allowed edit distance.
numThreads : int | None, default 4
    Number of worker threads for batch processing.

Returns
-------
list[list[str]]
    For each query, the list of matching sequences.
)pbdoc")

        .def("SearchAIRR",
             [](Trie& self,
                const std::string& query,
                int maxSubstitution,
                int maxInsertion,
                int maxDeletion,
                std::optional<int> maxEdits,
                const std::optional<std::string>& vGeneFilter,
                const std::optional<std::string>& jGeneFilter) {
                 return call_without_gil(&Trie::SearchAIRR,
                                         &self,
                                         query,
                                         maxSubstitution,
                                         maxInsertion,
                                         maxDeletion,
                                         maxEdits,
                                         vGeneFilter,
                                         jGeneFilter);
             },
             py::arg("query"),
             py::arg("maxSubstitution") = 0,
             py::arg("maxInsertion") = 0,
             py::arg("maxDeletion") = 0,
             py::arg("maxEdits") = std::nullopt,
             py::arg("vGeneFilter") = std::nullopt,
             py::arg("jGeneFilter") = std::nullopt,
             R"pbdoc(
Search and return AIRR-style results with optional V/J gene filtering.

Parameters
----------
query : str
    Query amino acid sequence.
maxSubstitution : int, default 0
    Maximum allowed substitutions.
maxInsertion : int, default 0
    Maximum allowed insertions.
maxDeletion : int, default 0
    Maximum allowed deletions.
maxEdits : int | None, default None
    Optional cap on total number of edits.
vGeneFilter : str | None, default None
    Restrict matches to this V gene.
jGeneFilter : str | None, default None
    Restrict matches to this J gene.

Returns
-------
list[AIRREntity]
    Matching AIRR-style records.
)pbdoc")

        .def("SearchForAll",
             [](Trie& self,
                const std::vector<std::string>& queries,
                int maxSubstitution,
                int maxInsertion,
                int maxDeletion,
                std::optional<int> maxEdits,
                std::optional<std::vector<std::string>> vGeneFilters,
                std::optional<std::vector<std::string>> jGeneFilters,
                std::optional<std::size_t> numThreads) {
                 return call_without_gil(&Trie::SearchForAll,
                                         &self,
                                         queries,
                                         maxSubstitution,
                                         maxInsertion,
                                         maxDeletion,
                                         maxEdits,
                                         vGeneFilters,
                                         jGeneFilters,
                                         numThreads);
             },
             py::arg("queries"),
             py::arg("maxSubstitution") = 0,
             py::arg("maxInsertion") = 0,
             py::arg("maxDeletion") = 0,
             py::arg("maxEdits") = std::nullopt,
             py::arg("vGeneFilters") = std::nullopt,
             py::arg("jGeneFilters") = std::nullopt,
             py::arg("numThreads") = 4,
             R"pbdoc(
Batch version of SearchAIRR.

Parameters
----------
queries : list[str]
    Query amino acid sequences.
maxSubstitution : int, default 0
    Maximum allowed substitutions.
maxInsertion : int, default 0
    Maximum allowed insertions.
maxDeletion : int, default 0
    Maximum allowed deletions.
maxEdits : int | None, default None
    Optional cap on total number of edits.
vGeneFilters : list[str] | None, default None
    Optional per-query V gene filters.
jGeneFilters : list[str] | None, default None
    Optional per-query J gene filters.
numThreads : int | None, default 4
    Number of worker threads.

Returns
-------
list[list[AIRREntity]]
    AIRR-style results for each query.
)pbdoc")

        .def("SearchWithMatrix",
             [](Trie& self,
                const std::string& query,
                float maxCost,
                const std::optional<std::string>& vGeneFilter,
                const std::optional<std::string>& jGeneFilter) {
                 return call_without_gil(&Trie::SearchWithMatrix,
                                         &self,
                                         query,
                                         maxCost,
                                         vGeneFilter,
                                         jGeneFilter);
             },
             py::arg("query"),
             py::arg("maxCost"),
             py::arg("vGeneFilter") = std::nullopt,
             py::arg("jGeneFilter") = std::nullopt,
             R"pbdoc(
Search using the loaded substitution matrix and a maximum total cost.

Parameters
----------
query : str
    Query amino acid sequence.
maxCost : float
    Maximum allowed alignment cost.
vGeneFilter : str | None, default None
    Restrict matches to this V gene.
jGeneFilter : str | None, default None
    Restrict matches to this J gene.

Returns
-------
list[AIRREntity]
    Matching AIRR-style records.
)pbdoc")

        .def("SearchForAllWithMatrix",
             [](Trie& self,
                const std::vector<std::string>& queries,
                float maxCost,
                std::optional<std::vector<std::string>> vGeneFilters,
                std::optional<std::vector<std::string>> jGeneFilters,
                std::optional<std::size_t> numThreads) {
                 return call_without_gil(&Trie::SearchForAllWithMatrix,
                                         &self,
                                         queries,
                                         maxCost,
                                         vGeneFilters,
                                         jGeneFilters,
                                         numThreads);
             },
             py::arg("queries"),
             py::arg("maxCost"),
             py::arg("vGeneFilters") = std::nullopt,
             py::arg("jGeneFilters") = std::nullopt,
             py::arg("numThreads") = 4,
             R"pbdoc(
Batch matrix-based search.

Parameters
----------
queries : list[str]
    Query amino acid sequences.
maxCost : float
    Maximum allowed alignment cost.
vGeneFilters : list[str] | None, default None
    Optional per-query V gene filters.
jGeneFilters : list[str] | None, default None
    Optional per-query J gene filters.
numThreads : int | None, default 4
    Number of worker threads.

Returns
-------
list[list[AIRREntity]]
    AIRR-style results for each query.
)pbdoc")

        .def("SearchIndices",
             [](Trie& self,
                const std::string& query,
                int maxSubstitution,
                int maxInsertion,
                int maxDeletion,
                std::optional<int> maxEdits,
                const std::optional<std::string>& vGeneFilter,
                const std::optional<std::string>& jGeneFilter) {
                 return call_without_gil(&Trie::SearchIndices,
                                         &self,
                                         query,
                                         maxSubstitution,
                                         maxInsertion,
                                         maxDeletion,
                                         maxEdits,
                                         vGeneFilter,
                                         jGeneFilter);
             },
             py::arg("query"),
             py::arg("maxSubstitution") = 0,
             py::arg("maxInsertion") = 0,
             py::arg("maxDeletion") = 0,
             py::arg("maxEdits") = std::nullopt,
             py::arg("vGeneFilter") = std::nullopt,
             py::arg("jGeneFilter") = std::nullopt,
             R"pbdoc(
Search and return index-based hits.

Parameters
----------
query : str
    Query amino acid sequence.
maxSubstitution : int, default 0
    Maximum allowed substitutions.
maxInsertion : int, default 0
    Maximum allowed insertions.
maxDeletion : int, default 0
    Maximum allowed deletions.
maxEdits : int | None, default None
    Optional cap on total edits.
vGeneFilter : str | None, default None
    Restrict matches to this V gene.
jGeneFilter : str | None, default None
    Restrict matches to this J gene.

Returns
-------
list[tuple[int, int]]
    Pairs (target_index, distance).
)pbdoc")

        .def("SearchIndicesForAll",
             [](Trie& self,
                const std::vector<std::string>& queries,
                int maxSubstitution,
                int maxInsertion,
                int maxDeletion,
                std::optional<int> maxEdits,
                std::optional<std::vector<std::string>> vGeneFilters,
                std::optional<std::vector<std::string>> jGeneFilters,
                std::optional<std::size_t> numThreads) {
                 return call_without_gil(&Trie::SearchIndicesForAll,
                                         &self,
                                         queries,
                                         maxSubstitution,
                                         maxInsertion,
                                         maxDeletion,
                                         maxEdits,
                                         vGeneFilters,
                                         jGeneFilters,
                                         numThreads);
             },
             py::arg("queries"),
             py::arg("maxSubstitution") = 0,
             py::arg("maxInsertion") = 0,
             py::arg("maxDeletion") = 0,
             py::arg("maxEdits") = std::nullopt,
             py::arg("vGeneFilters") = std::nullopt,
             py::arg("jGeneFilters") = std::nullopt,
             py::arg("numThreads") = 4,
             R"pbdoc(
Batch version of SearchIndices.

Parameters
----------
queries : list[str]
    Query amino acid sequences.
maxSubstitution : int, default 0
    Maximum allowed substitutions.
maxInsertion : int, default 0
    Maximum allowed insertions.
maxDeletion : int, default 0
    Maximum allowed deletions.
maxEdits : int | None, default None
    Optional cap on total edits.
vGeneFilters : list[str] | None, default None
    Optional per-query V gene filters.
jGeneFilters : list[str] | None, default None
    Optional per-query J gene filters.
numThreads : int | None, default 4
    Number of worker threads.

Returns
-------
list[list[tuple[int, int]]]
    Index-based hits for each query.
)pbdoc")

        .def("SearchIndicesWithMatrix",
             [](Trie& self,
                const std::string& query,
                float maxCost,
                const std::optional<std::string>& vGeneFilter,
                const std::optional<std::string>& jGeneFilter) {
                 return call_without_gil(&Trie::SearchIndicesWithMatrix,
                                         &self,
                                         query,
                                         maxCost,
                                         vGeneFilter,
                                         jGeneFilter);
             },
             py::arg("query"),
             py::arg("maxCost"),
             py::arg("vGeneFilter") = std::nullopt,
             py::arg("jGeneFilter") = std::nullopt,
             R"pbdoc(
Matrix-based search returning index-based hits.

Parameters
----------
query : str
    Query amino acid sequence.
maxCost : float
    Maximum allowed alignment cost.
vGeneFilter : str | None, default None
    Restrict matches to this V gene.
jGeneFilter : str | None, default None
    Restrict matches to this J gene.

Returns
-------
list[tuple[int, float]]
    Pairs (target_index, cost).
)pbdoc")

        .def("SearchIndicesForAllWithMatrix",
             [](Trie& self,
                const std::vector<std::string>& queries,
                float maxCost,
                std::optional<std::vector<std::string>> vGeneFilters,
                std::optional<std::vector<std::string>> jGeneFilters,
                std::optional<std::size_t> numThreads) {
                 return call_without_gil(&Trie::SearchIndicesForAllWithMatrix,
                                         &self,
                                         queries,
                                         maxCost,
                                         vGeneFilters,
                                         jGeneFilters,
                                         numThreads);
             },
             py::arg("queries"),
             py::arg("maxCost"),
             py::arg("vGeneFilters") = std::nullopt,
             py::arg("jGeneFilters") = std::nullopt,
             py::arg("numThreads") = 4,
             R"pbdoc(
Batch matrix-based search returning index-based hits.

Parameters
----------
queries : list[str]
    Query amino acid sequences.
maxCost : float
    Maximum allowed alignment cost.
vGeneFilters : list[str] | None, default None
    Optional per-query V gene filters.
jGeneFilters : list[str] | None, default None
    Optional per-query J gene filters.
numThreads : int | None, default 4
    Number of worker threads.

Returns
-------
list[list[tuple[int, float]]]
    Index-based hits for each query.
)pbdoc")

        .def("SearchGroupIdsForAll",
             [](Trie& self,
                const std::vector<std::string>& queries,
                int maxSubstitution,
                int maxInsertion,
                int maxDeletion,
                std::optional<int> maxEdits,
                std::optional<std::vector<std::string>> vGeneFilters,
                std::optional<std::vector<std::string>> jGeneFilters,
                bool unique,
                std::optional<std::size_t> numThreads) {
                 return call_without_gil(&Trie::SearchGroupIdsForAll,
                                         &self,
                                         queries,
                                         maxSubstitution,
                                         maxInsertion,
                                         maxDeletion,
                                         maxEdits,
                                         vGeneFilters,
                                         jGeneFilters,
                                         unique,
                                         numThreads);
             },
             py::arg("queries"),
             py::arg("maxSubstitution") = 0,
             py::arg("maxInsertion") = 0,
             py::arg("maxDeletion") = 0,
             py::arg("maxEdits") = std::nullopt,
             py::arg("vGeneFilters") = std::nullopt,
             py::arg("jGeneFilters") = std::nullopt,
             py::arg("unique") = true,
             py::arg("numThreads") = 4,
             R"pbdoc(
For each query, return matching group identifiers.

Parameters
----------
queries : list[str]
    Query amino acid sequences.
maxSubstitution : int, default 0
    Maximum allowed substitutions.
maxInsertion : int, default 0
    Maximum allowed insertions.
maxDeletion : int, default 0
    Maximum allowed deletions.
maxEdits : int | None, default None
    Optional cap on total edits.
vGeneFilters : list[str] | None, default None
    Optional per-query V gene filters.
jGeneFilters : list[str] | None, default None
    Optional per-query J gene filters.
unique : bool, default True
    Whether to deduplicate group ids in the result.
numThreads : int | None, default 4
    Number of worker threads.

Returns
-------
object
    Group id results for each query.
)pbdoc")

        .def("ClusterUsage",
             [](Trie& self,
                const std::vector<std::string>& cluster,
                int maxSubstitution,
                int maxInsertion,
                int maxDeletion,
                std::optional<int> maxEdits,
                std::optional<std::vector<std::string>> vGeneFilters,
                std::optional<std::vector<std::string>> jGeneFilters,
                std::optional<std::size_t> numThreads) {
                 return call_without_gil(&Trie::ClusterUsage,
                                         &self,
                                         cluster,
                                         maxSubstitution,
                                         maxInsertion,
                                         maxDeletion,
                                         maxEdits,
                                         vGeneFilters,
                                         jGeneFilters,
                                         numThreads);
             },
             py::arg("cluster"),
             py::arg("maxSubstitution") = 0,
             py::arg("maxInsertion") = 0,
             py::arg("maxDeletion") = 0,
             py::arg("maxEdits") = std::nullopt,
             py::arg("vGeneFilters") = std::nullopt,
             py::arg("jGeneFilters") = std::nullopt,
             py::arg("numThreads") = 4,
             R"pbdoc(
Estimate usage of a sequence cluster with edit-distance-based matching.

Parameters
----------
cluster : list[str]
    Cluster sequences to evaluate.
maxSubstitution : int, default 0
    Maximum allowed substitutions.
maxInsertion : int, default 0
    Maximum allowed insertions.
maxDeletion : int, default 0
    Maximum allowed deletions.
maxEdits : int | None, default None
    Optional cap on total edits.
vGeneFilters : list[str] | None, default None
    Optional per-sequence V gene filters.
jGeneFilters : list[str] | None, default None
    Optional per-sequence J gene filters.
numThreads : int | None, default 4
    Number of worker threads.

Returns
-------
object
    Cluster usage summary.
)pbdoc")

        .def("ClusterUsageWithMatrix",
             [](Trie& self,
                const std::vector<std::string>& cluster,
                float maxCost,
                std::optional<std::vector<std::string>> vGeneFilters,
                std::optional<std::vector<std::string>> jGeneFilters,
                std::optional<std::size_t> numThreads) {
                 return call_without_gil(&Trie::ClusterUsageWithMatrix,
                                         &self,
                                         cluster,
                                         maxCost,
                                         vGeneFilters,
                                         jGeneFilters,
                                         numThreads);
             },
             py::arg("cluster"),
             py::arg("maxCost"),
             py::arg("vGeneFilters") = std::nullopt,
             py::arg("jGeneFilters") = std::nullopt,
             py::arg("numThreads") = 4,
             R"pbdoc(
Estimate usage of a sequence cluster with matrix-based matching.

Parameters
----------
cluster : list[str]
    Cluster sequences to evaluate.
maxCost : float
    Maximum allowed alignment cost.
vGeneFilters : list[str] | None, default None
    Optional per-sequence V gene filters.
jGeneFilters : list[str] | None, default None
    Optional per-sequence J gene filters.
numThreads : int | None, default 4
    Number of worker threads.

Returns
-------
object
    Cluster usage summary.
)pbdoc")

        .def("SearchAny",
             [](Trie& self, const std::string& query, int maxEdits) {
                 return call_without_gil(&Trie::SearchAny, &self, query, maxEdits);
             },
             py::arg("query"),
             py::arg("maxEdits"),
             R"pbdoc(
Check whether at least one match exists within the given edit distance.

Parameters
----------
query : str
    Query amino acid sequence.
maxEdits : int
    Maximum allowed edit distance.

Returns
-------
bool
    True if at least one match exists, otherwise False.
)pbdoc")

        .def("AlignQueryToTarget",
             [](Trie& self,
                const std::string& query,
                const std::string& target,
                std::optional<int> maxSubstitution,
                std::optional<int> maxInsertion,
                std::optional<int> maxDeletion,
                std::optional<int> maxEdits) {
                 return call_without_gil(&Trie::AlignQueryToTarget,
                                         &self,
                                         query,
                                         target,
                                         maxSubstitution,
                                         maxInsertion,
                                         maxDeletion,
                                         maxEdits);
             },
             py::arg("query"),
             py::arg("target"),
             py::arg("maxSubstitution") = std::nullopt,
             py::arg("maxInsertion") = std::nullopt,
             py::arg("maxDeletion") = std::nullopt,
             py::arg("maxEdits") = std::nullopt,
             R"pbdoc(
Build a detailed alignment between a query and an explicit target sequence.

Parameters
----------
query : str
    Query amino acid sequence.
target : str
    Target amino acid sequence.
maxSubstitution : int | None, default None
    Optional substitution limit.
maxInsertion : int | None, default None
    Optional insertion limit.
maxDeletion : int | None, default None
    Optional deletion limit.
maxEdits : int | None, default None
    Optional total edit limit.

Returns
-------
AlignmentResult
    Detailed alignment information.
)pbdoc")

        .def("AlignIndexHit",
             [](Trie& self,
                const std::string& query,
                size_t targetIndex,
                std::optional<int> maxSubstitution,
                std::optional<int> maxInsertion,
                std::optional<int> maxDeletion,
                std::optional<int> maxEdits) {
                 return call_without_gil(&Trie::AlignIndexHit,
                                         &self,
                                         query,
                                         targetIndex,
                                         maxSubstitution,
                                         maxInsertion,
                                         maxDeletion,
                                         maxEdits);
             },
             py::arg("query"),
             py::arg("targetIndex"),
             py::arg("maxSubstitution") = std::nullopt,
             py::arg("maxInsertion") = std::nullopt,
             py::arg("maxDeletion") = std::nullopt,
             py::arg("maxEdits") = std::nullopt,
             R"pbdoc(
Build a detailed alignment between a query and a target stored in the trie by index.

Parameters
----------
query : str
    Query amino acid sequence.
targetIndex : int
    Index of the stored target sequence.
maxSubstitution : int | None, default None
    Optional substitution limit.
maxInsertion : int | None, default None
    Optional insertion limit.
maxDeletion : int | None, default None
    Optional deletion limit.
maxEdits : int | None, default None
    Optional total edit limit.

Returns
-------
AlignmentResult
    Detailed alignment information.
)pbdoc")

        .def("AlignIndexHits",
             [](Trie& self,
                const std::string& query,
                const std::vector<std::pair<size_t, int>>& hits,
                std::optional<int> maxSubstitution,
                std::optional<int> maxInsertion,
                std::optional<int> maxDeletion,
                std::optional<int> maxEdits,
                std::optional<std::size_t> numThreads) {
                 return call_without_gil(&Trie::AlignIndexHits,
                                         &self,
                                         query,
                                         hits,
                                         maxSubstitution,
                                         maxInsertion,
                                         maxDeletion,
                                         maxEdits,
                                         numThreads);
             },
             py::arg("query"),
             py::arg("hits"),
             py::arg("maxSubstitution") = std::nullopt,
             py::arg("maxInsertion") = std::nullopt,
             py::arg("maxDeletion") = std::nullopt,
             py::arg("maxEdits") = std::nullopt,
             py::arg("numThreads") = 4,
             R"pbdoc(
Build detailed alignments for multiple index-based edit-distance hits.

Parameters
----------
query : str
    Query amino acid sequence.
hits : list[tuple[int, int]]
    Pairs (target_index, distance), usually returned by SearchIndices.
maxSubstitution : int | None, default None
    Optional substitution limit.
maxInsertion : int | None, default None
    Optional insertion limit.
maxDeletion : int | None, default None
    Optional deletion limit.
maxEdits : int | None, default None
    Optional total edit limit.
numThreads : int | None, default 4
    Number of worker threads.

Returns
-------
list[AlignmentResult]
    Detailed alignments for the provided hits.
)pbdoc")

        .def("AlignQueryToTargetWithMatrix",
             [](Trie& self,
                const std::string& query,
                const std::string& target,
                std::optional<float> maxCost) {
                 return call_without_gil(&Trie::AlignQueryToTargetWithMatrix,
                                         &self,
                                         query,
                                         target,
                                         maxCost);
             },
             py::arg("query"),
             py::arg("target"),
             py::arg("maxCost") = std::nullopt,
             R"pbdoc(
Build a detailed alignment using the loaded substitution matrix.

Parameters
----------
query : str
    Query amino acid sequence.
target : str
    Target amino acid sequence.
maxCost : float | None, default None
    Optional maximum allowed alignment cost.

Returns
-------
AlignmentResult
    Detailed alignment information.
)pbdoc")

        .def("AlignIndexHitWithMatrix",
             [](Trie& self,
                const std::string& query,
                size_t targetIndex,
                std::optional<float> maxCost) {
                 return call_without_gil(&Trie::AlignIndexHitWithMatrix,
                                         &self,
                                         query,
                                         targetIndex,
                                         maxCost);
             },
             py::arg("query"),
             py::arg("targetIndex"),
             py::arg("maxCost") = std::nullopt,
             R"pbdoc(
Build a detailed matrix-based alignment against a stored target by index.

Parameters
----------
query : str
    Query amino acid sequence.
targetIndex : int
    Index of the stored target sequence.
maxCost : float | None, default None
    Optional maximum allowed alignment cost.

Returns
-------
AlignmentResult
    Detailed alignment information.
)pbdoc")

        .def("AlignIndexHitsWithMatrix",
             [](Trie& self,
                const std::string& query,
                const std::vector<std::pair<size_t, float>>& hits,
                std::optional<float> maxCost,
                std::optional<std::size_t> numThreads) {
                 return call_without_gil(&Trie::AlignIndexHitsWithMatrix,
                                         &self,
                                         query,
                                         hits,
                                         maxCost,
                                         numThreads);
             },
             py::arg("query"),
             py::arg("hits"),
             py::arg("maxCost") = std::nullopt,
             py::arg("numThreads") = 4,
             R"pbdoc(
Build detailed alignments for multiple index-based matrix hits.

Parameters
----------
query : str
    Query amino acid sequence.
hits : list[tuple[int, float]]
    Pairs (target_index, cost), usually returned by SearchIndicesWithMatrix.
maxCost : float | None, default None
    Optional maximum allowed alignment cost.
numThreads : int | None, default 4
    Number of worker threads.

Returns
-------
list[AlignmentResult]
    Detailed alignments for the provided hits.
)pbdoc");
}
