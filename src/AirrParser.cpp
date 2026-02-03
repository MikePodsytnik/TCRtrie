#include "AirrParser.h"

#include <fstream>
#include <iostream>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <cctype>

namespace {
    inline bool IsWhitespace(char c) {
        return std::isspace(static_cast<unsigned char>(c));
    }

    size_t FindNextWhitespace(const std::string& str, size_t pos) {
        while (pos < str.size() && !IsWhitespace(str[pos])) {
            ++pos;
        }
        return pos;
    }

    size_t SkipWhitespace(const std::string& str, size_t pos) {
        while (pos < str.size() && IsWhitespace(str[pos])) {
            ++pos;
        }
        return pos;
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

    std::unordered_map<std::string, int> colIdx;
    {
        int idx = 0;
        size_t pos = 0;

        while (pos < header.size()) {
            pos = SkipWhitespace(header, pos);
            if (pos >= header.size()) break;

            size_t end = FindNextWhitespace(header, pos);

            std::string colName = header.substr(pos, end - pos);
            if (!colName.empty()) {
                colIdx[colName] = idx++;
            }

            pos = end;
        }
    }

    auto itJ = colIdx.find("junction_aa");
    if (itJ == colIdx.end()) {
        std::cerr << "[Error] No column junction_aa.\n";
        return entries;
    }
    int junctionCol = itJ->second;
    int vCol = colIdx.count("v_call") ? colIdx["v_call"] : -1;
    int jCol = colIdx.count("j_call") ? colIdx["j_call"] : -1;
    int maxCol = std::max(junctionCol, std::max(vCol, jCol));

    std::string line;
    while (std::getline(file, line)) {
        AIRREntity ent;
        int col = 0;
        size_t start = 0;

        while (col <= maxCol && start < line.size()) {
            start = SkipWhitespace(line, start);
            if (start >= line.size()) break;

            size_t end = FindNextWhitespace(line, start);

            std::string_view fv{ line.data() + start, end - start };

            if (col == junctionCol) {
                ent.junctionAA.assign(fv);
            }
            else if (col == vCol) {
                ent.vGene.assign(fv);
            }
            else if (col == jCol) {
                ent.jGene.assign(fv);
            }

            start = end;
            ++col;
        }

        if (!ent.junctionAA.empty()) {
            entries.push_back(std::move(ent));
        }
    }

    return entries;
}
