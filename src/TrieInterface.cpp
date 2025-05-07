#include "TrieInterface.h"
#include "Trie.h"

#include <fstream>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

static void WriteResults(const std::string& outPath, const std::unordered_map<std::string, std::vector<AIRREntity>>& results) {
    std::ofstream outFile(outPath);
    if (!outFile.is_open()) {
        std::cerr << "Error: Unable to write to " << outPath << std::endl;
        return;
    }

    bool hasVGene = false, hasJGene = false;
    for (const auto& [_, matches] : results) {
        for (const auto& m : matches) {
            if (!m.vGene.empty()) hasVGene = true;
            if (!m.jGene.empty()) hasJGene = true;
            if (hasVGene && hasJGene) break;
        }
        if (hasVGene && hasJGene) break;
    }

    outFile << "query\tmatch";
    if (hasVGene) outFile << "\tv_gene";
    if (hasJGene) outFile << "\tj_gene";
    outFile << '\n';

    for (const auto& [query, matches] : results) {
        for (const auto& match : matches) {
            outFile << query << '\t' << match.junctionAA;
            if (hasVGene) outFile << '\t' << match.vGene;
            if (hasJGene) outFile << '\t' << match.jGene;
            outFile << '\n';
        }
    }
}

static std::vector<std::string> LoadQueriesFromFile(const std::string& path) {
    std::ifstream file(path);
    std::vector<std::string> queries;
    std::string line;

    std::getline(file, line);

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string field;
        std::getline(ss, field, '\t');
        if (!field.empty()) {
            queries.push_back(field);
        }
    }
    return queries;
}

void RunSearch(const SearchConfig& config) {
    Trie trie(config.inputPath);

    if (!config.matrixPath.empty()) {
        trie.LoadSubstitutionMatrix(config.matrixPath);
    }

    fs::create_directories(config.outputPath);
    std::string outFilePath = config.outputPath + "/results.tsv";

    if (!config.query.empty()) {
        std::vector<AIRREntity> results;

        if (!config.matrixPath.empty()) {
            results = trie.SearchWithScore(config.query, config.scoreRadius);
        } else {
            results = trie.SearchAIRR(config.query, config.nEdits);
        }

        std::unordered_map<std::string, std::vector<AIRREntity>> wrapped{
                {config.query, results}
        };
        WriteResults(outFilePath, wrapped);
    }
    else if (!config.inputQueries.empty()) {
        std::vector<std::string> queries = LoadQueriesFromFile(config.inputQueries);

        std::unordered_map<std::string, std::vector<AIRREntity>> resultMap;

        if (!config.matrixPath.empty()) {
            resultMap = trie.SearchForAllWithScore(queries, config.scoreRadius);
        } else {
            resultMap = trie.SearchForAll(queries, config.nEdits);
        }

        WriteResults(outFilePath, resultMap);
    }
    else {
        std::cerr << "Error: No query provided.\n";
    }

    std::cout << "SearchAIRR complete. Results saved to: " << outFilePath << std::endl;
}
