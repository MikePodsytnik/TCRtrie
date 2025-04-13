#include "AirrParser.h"
#include <fstream>
#include <sstream>
#include <iostream>

std::vector<AIRREntity> ParseAIRR(const std::string& filepath) {
    std::vector<AIRREntity> entries;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "[Error] Failed to open AIRR file: " << filepath << '\n';
        return entries;
    }

    std::string headerLine;
    if (!std::getline(file, headerLine)) {
        std::cerr << "[Error] AIRR file is empty.\n";
        return entries;
    }

    std::unordered_map<std::string, int> columnIndices;
    std::istringstream headerStream(headerLine);
    std::string column;
    int idx = 0;
    while (std::getline(headerStream, column, '\t')) {
        columnIndices[column] = idx++;
    }

    if (columnIndices.find("junction_aa") == columnIndices.end()) {
        std::cerr << "[Error] Required column junction_aa not found.\n";
        return entries;
    }

    int junctionCol = columnIndices["junction_aa"];
    int vCol = columnIndices.count("v_call") ? columnIndices["v_call"] : -1;
    int jCol = columnIndices.count("j_call") ? columnIndices["j_call"] : -1;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::vector<std::string> fields;
        std::string field;

        while (std::getline(lineStream, field, '\t')) {
            fields.push_back(std::move(field));
        }

        if (fields.size() <= junctionCol) continue;

        AIRREntity entry;
        entry.junctionAA = fields[junctionCol];
        if (entry.junctionAA.empty()) continue;

        if (vCol != -1 && vCol < fields.size()) entry.vGene = fields[vCol];
        if (jCol != -1 && jCol < fields.size()) entry.jGene = fields[jCol];

        entries.push_back(std::move(entry));
    }

    return entries;
}
