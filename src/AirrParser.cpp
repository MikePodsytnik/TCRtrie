#include "AirrParser.h"

#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <charconv>

namespace {

inline void TrimCR(std::string& s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
}

inline std::vector<std::string_view> SplitTSVViews(const std::string& line) {
    std::vector<std::string_view> fields;
    fields.reserve(32);

    size_t start = 0;
    while (start <= line.size()) {
        size_t end = line.find('\t', start);
        if (end == std::string::npos) end = line.size();
        fields.emplace_back(line.data() + start, end - start);
        start = end + 1;
        if (end == line.size()) break;
    }
    return fields;
}

inline int ParseIntOrDefault(std::string_view sv, int def) {
    if (sv.empty()) return def;

    int value = def;
    auto first = sv.data();
    auto last = sv.data() + sv.size();

    while (first < last && (*first == ' ' || *first == '\t')) ++first;
    while (last > first && (*(last - 1) == ' ' || *(last - 1) == '\t')) --last;

    if (first >= last) return def;

    std::from_chars_result r = std::from_chars(first, last, value);
    if (r.ec != std::errc()) return def;
    return value;
}

}

std::vector<AIRREntity> ParseAIRR(const std::string& filepath) {
    std::vector<AIRREntity> entries;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "[Error] Failed to open " << filepath << '\n';
        return entries;
    }

    std::string header;
    if (!std::getline(file, header)) {
        std::cerr << "[Error] Empty file.\n";
        return entries;
    }
    TrimCR(header);

    std::unordered_map<std::string, int> colIdx;
    {
        auto fields = SplitTSVViews(header);
        colIdx.reserve(fields.size());

        for (int i = 0; i < static_cast<int>(fields.size()); ++i) {
            std::string name(fields[i]);
            if (!name.empty()) colIdx.emplace(std::move(name), i);
        }
    }

    auto itJ = colIdx.find("junction_aa");
    if (itJ == colIdx.end()) {
        std::cerr << "[Error] No column junction_aa.\n";
        return entries;
    }

    const int junctionCol = itJ->second;
    const int vCol = colIdx.count("v_call") ? colIdx["v_call"] : -1;
    const int jCol = colIdx.count("j_call") ? colIdx["j_call"] : -1;
    const int gidCol = colIdx.count("__group_id") ? colIdx["__group_id"] : -1;

    int maxCol = junctionCol;
    if (vCol > maxCol) maxCol = vCol;
    if (jCol > maxCol) maxCol = jCol;
    if (gidCol > maxCol) maxCol = gidCol;

    std::string line;
    while (std::getline(file, line)) {
        TrimCR(line);
        if (line.empty()) continue;

        auto fields = SplitTSVViews(line);
        if (junctionCol >= static_cast<int>(fields.size())) continue;

        AIRREntity ent;

        {
            auto sv = fields[junctionCol];
            if (sv.empty()) continue;
            ent.junctionAA.assign(sv);
        }

        if (vCol >= 0 && vCol < static_cast<int>(fields.size())) {
            auto sv = fields[vCol];
            if (!sv.empty()) ent.vGene.assign(sv);
        }

        if (jCol >= 0 && jCol < static_cast<int>(fields.size())) {
            auto sv = fields[jCol];
            if (!sv.empty()) ent.jGene.assign(sv);
        }

        if (gidCol >= 0 && gidCol < static_cast<int>(fields.size())) {
            ent.groupId = ParseIntOrDefault(fields[gidCol], -1);
        } else {
            ent.groupId = -1;
        }

        entries.push_back(std::move(ent));
    }

    return entries;
}