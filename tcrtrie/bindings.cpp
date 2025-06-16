#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
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
            .def("__repr__", [](const AIRREntity &e) {
                return "<AIRREntity junctionAA='" + e.junctionAA + "' vGene='" + e.vGene +
                       "' jGene='" + e.jGene + "' distance=" + std::to_string(e.distance) + ">";
            });

    py::class_<Trie>(m, "Trie")
            .def(py::init<>())
            .def(py::init<const std::vector<std::string>&>(), py::arg("sequences"))
            .def(py::init<const std::string&>(), py::arg("dataPath"))
            .def(py::init<const Trie&>(), py::arg("other"))

            .def("Search",
                 [](Trie &self, const std::string &query, int maxEdits) {
                     using Method = std::vector<std::string> (Trie::*)(const std::string&, int);
                     return call_without_gil(static_cast<Method>(&Trie::Search), &self, query, maxEdits);
                 },
                 py::arg("query"), py::arg("maxEdits"))

            .def("Search",
                 [](Trie &self, const std::vector<std::string> &queries, int maxEdits) {
                     using Method = std::unordered_map<std::string, std::vector<std::string>> (Trie::*)(const std::vector<std::string>&, int);
                     return call_without_gil(static_cast<Method>(&Trie::Search), &self, queries, maxEdits);
                 },
                 py::arg("queries"), py::arg("maxEdits"))

            .def("SearchAIRR",
                 [](Trie &self,
                    const std::string &query,
                    int maxSubstitution,
                    int maxInsertion,
                    int maxDeletion,
                    const std::optional<std::string> &vGeneFilter,
                    const std::optional<std::string> &jGeneFilter) {
                     return call_without_gil(&Trie::SearchAIRR,
                                             &self,
                                             query,
                                             maxSubstitution,
                                             maxInsertion,
                                             maxDeletion,
                                             vGeneFilter,
                                             jGeneFilter);
                 },
                 py::arg("query"),
                 py::arg("maxSubstitution"),
                 py::arg("maxInsertion"),
                 py::arg("maxDeletion"),
                 py::arg("vGeneFilter") = std::nullopt,
                 py::arg("jGeneFilter") = std::nullopt)

            .def("SearchWithMatrix",
                 [](Trie &self,
                    const std::string &query,
                    float maxCost,
                    const std::optional<std::string> &vGeneFilter,
                    const std::optional<std::string> &jGeneFilter) {
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

            .def("SearchAny",
                 [](Trie &self, const std::string &query, int maxEdits) {
                     return call_without_gil(&Trie::SearchAny, &self, query, maxEdits);
                 },
                 py::arg("query"), py::arg("maxEdits"))

            .def("SearchForAll",
                 [](Trie &self,
                    const std::vector<std::string> &queries,
                    int maxSubstitution,
                    int maxInsertion,
                    int maxDeletion,
                    const std::optional<std::string> &vGeneFilter,
                    const std::optional<std::string> &jGeneFilter) {
                     return call_without_gil(&Trie::SearchForAll,
                                             &self,
                                             queries,
                                             maxSubstitution,
                                             maxInsertion,
                                             maxDeletion,
                                             vGeneFilter,
                                             jGeneFilter);
                 },
                 py::arg("queries"),
                 py::arg("maxSubstitution"),
                 py::arg("maxInsertion"),
                 py::arg("maxDeletion"),
                 py::arg("vGeneFilter") = std::nullopt,
                 py::arg("jGeneFilter") = std::nullopt)

            .def("SearchForAllWithMatrix",
                 [](Trie &self,
                    const std::vector<std::string> &queries,
                    float maxCost,
                    const std::optional<std::string> &vGeneFilter,
                    const std::optional<std::string> &jGeneFilter) {
                     return call_without_gil(&Trie::SearchForAllWithMatrix,
                                             &self,
                                             queries,
                                             maxCost,
                                             vGeneFilter,
                                             jGeneFilter);
                 },
                 py::arg("queries"),
                 py::arg("maxCost"),
                 py::arg("vGeneFilter") = std::nullopt,
                 py::arg("jGeneFilter") = std::nullopt)

            .def("LoadSubstitutionMatrix", &Trie::LoadSubstitutionMatrix, py::arg("matrixPath"))

            .def("SetMaxQueryLength", &Trie::SetMaxQueryLength, py::arg("newMaxQueryLength"));
}