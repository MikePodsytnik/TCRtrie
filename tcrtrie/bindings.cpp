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

    py::class_<AIRREntity>(m, "AIRREntity")
        .def_readonly("junctionAA", &AIRREntity::junctionAA)
        .def_readonly("vGene", &AIRREntity::vGene)
        .def_readonly("jGene", &AIRREntity::jGene)
        .def_readonly("distance", &AIRREntity::distance)
        .def("__repr__", [](const AIRREntity& e) {
            return "<AIRREntity junctionAA='" + e.junctionAA + "' vGene='" + e.vGene +
                   "' jGene='" + e.jGene + "' distance=" + std::to_string(e.distance) + ">";
        });

    py::enum_<Trie::AlignmentOpType>(m, "AlignmentOpType")
        .value("Match", Trie::AlignmentOpType::Match)
        .value("Substitution", Trie::AlignmentOpType::Substitution)
        .value("Insertion", Trie::AlignmentOpType::Insertion)
        .value("Deletion", Trie::AlignmentOpType::Deletion)
        .export_values();

    py::class_<Trie::AlignmentOp>(m, "AlignmentOp")
        .def_readonly("type", &Trie::AlignmentOp::type)
        .def_readonly("queryPos", &Trie::AlignmentOp::queryPos)
        .def_readonly("queryChar", &Trie::AlignmentOp::queryChar)
        .def_readonly("targetChar", &Trie::AlignmentOp::targetChar)
        .def("__repr__", [](const Trie::AlignmentOp& op) {
            return "<AlignmentOp type=" + py::str(py::cast(op.type)).cast<std::string>() +
                   " queryPos=" + std::to_string(op.queryPos) +
                   " queryChar='" + std::string(1, op.queryChar) +
                   "' targetChar='" + std::string(1, op.targetChar) + "'>";
        });

    py::class_<Trie::AlignmentResult>(m, "AlignmentResult")
        .def_readonly("queryAligned", &Trie::AlignmentResult::queryAligned)
        .def_readonly("targetAligned", &Trie::AlignmentResult::targetAligned)
        .def_readonly("substitutions", &Trie::AlignmentResult::substitutions)
        .def_readonly("insertions", &Trie::AlignmentResult::insertions)
        .def_readonly("deletions", &Trie::AlignmentResult::deletions)
        .def_readonly("distance", &Trie::AlignmentResult::distance)
        .def_readonly("ops", &Trie::AlignmentResult::ops)
        .def("__repr__", [](const Trie::AlignmentResult& r) {
            return "<AlignmentResult distance=" + std::to_string(r.distance) +
                   " substitutions=" + std::to_string(r.substitutions) +
                   " insertions=" + std::to_string(r.insertions) +
                   " deletions=" + std::to_string(r.deletions) +
                   " queryAligned='" + r.queryAligned +
                   "' targetAligned='" + r.targetAligned + "'>";
        });

    py::class_<Trie>(m, "Trie")
        .def(py::init<>())
        .def(py::init<const std::string&>(), py::arg("dataPath"))
        .def(py::init<const Trie&>(), py::arg("other"))
        .def(py::init<const std::vector<std::string>&,
                      const std::vector<std::string>&,
                      const std::vector<std::string>&>(),
             py::arg("sequences"),
             py::arg("vGenes"),
             py::arg("jGenes"))

        .def("LoadSubstitutionMatrix",
             &Trie::LoadSubstitutionMatrix,
             py::arg("matrixPath"),
             py::arg("delimiter") = "",
             py::arg("gapFactor") = 1.5f)
        .def("PrintMatrix", &Trie::PrintMatrix)

        .def("Search",
             [](Trie& self, const std::string& query, int maxEdits) {
                 using Method = std::vector<std::string> (Trie::*)(const std::string&, int);
                 return call_without_gil(static_cast<Method>(&Trie::Search), &self, query, maxEdits);
             },
             py::arg("query"), py::arg("maxEdits"))

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
             py::arg("numThreads") = 4)

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
             py::arg("jGeneFilter") = std::nullopt)

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
             py::arg("numThreads") = 4)

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
             py::arg("jGeneFilter") = std::nullopt)

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
             py::arg("numThreads") = 4)

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
             py::arg("jGeneFilter") = std::nullopt)

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
             py::arg("numThreads") = 4)

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
             py::arg("jGeneFilter") = std::nullopt)

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
             py::arg("numThreads") = 4)

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
             py::arg("numThreads") = 4)

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
             py::arg("numThreads") = 4)

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
             py::arg("numThreads") = 4)

        .def("SearchAny",
             [](Trie& self, const std::string& query, int maxEdits) {
                 return call_without_gil(&Trie::SearchAny, &self, query, maxEdits);
             },
             py::arg("query"), py::arg("maxEdits"))

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
             py::arg("maxEdits") = std::nullopt)

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
             py::arg("maxEdits") = std::nullopt)

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
             py::arg("numThreads") = 4)

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
             py::arg("maxCost") = std::nullopt)

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
             py::arg("maxCost") = std::nullopt)

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
             py::arg("numThreads") = 4);
}
