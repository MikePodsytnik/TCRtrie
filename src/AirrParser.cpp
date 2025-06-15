#include "AirrParser.h"
#include <fstream>
#include <iostream>
#include <string_view>

std::vector<AIRREntity> ParseAIRR(const std::string& filepath) {
    std::vector<AIRREntity> entries;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "[Error] Не удалось открыть " << filepath << '\n';
        return entries;
    }

    std::string header;
    if (!std::getline(file, header)) {
        std::cerr << "[Error] Пустой файл.\n";
        return entries;
    }
    std::unordered_map<std::string, int> colIdx;
    {
        int idx = 0, pos = 0;
        while (pos <= (int)header.size()) {
            int tab = header.find('\t', pos);
            if (tab == std::string::npos) tab = header.size();
            colIdx[ header.substr(pos, tab - pos) ] = idx++;
            pos = tab + 1;
        }
    }
    auto itJ = colIdx.find("junction_aa");
    if (itJ == colIdx.end()) {
        std::cerr << "[Error] Нет колонки junction_aa.\n";
        return entries;
    }
    int junctionCol = itJ->second;
    int vCol = colIdx.count("v_call") ? colIdx["v_call"] : -1;
    int jCol = colIdx.count("j_call") ? colIdx["j_call"] : -1;
    int maxCol = std::max({junctionCol, vCol, jCol});

    std::string line;
    while (std::getline(file, line)) {
        AIRREntity ent;
        int col = 0;
        size_t start = 0;

        while (col <= maxCol && start <= line.size()) {
            size_t end = line.find('\t', start);
            if (end == std::string::npos) end = line.size();

            std::string_view fv{ line.data() + start, end - start };
            if (col == junctionCol) {
                ent.junctionAA.assign(fv);
                if (ent.junctionAA.empty()) break;
            }
            else if (col == vCol) {
                ent.vGene.assign(fv);
            }
            else if (col == jCol) {
                ent.jGene.assign(fv);
            }

            start = end + 1;
            ++col;
        }

        if (!ent.junctionAA.empty()) {
            entries.push_back(std::move(ent));
        }
    }

    return entries;
}
