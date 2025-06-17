#include "TrieInterface.h"
#include "Trie.h"

#include <fstream>
#include <filesystem>
#include <sstream>
#include <iostream>
#include <algorithm>

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

    outFile << "query\tmatch\tdist";
    if (hasVGene) outFile << "\tv_gene";
    if (hasJGene) outFile << "\tj_gene";
    outFile << '\n';

    for (const auto& [query, matches] : results) {
        for (const auto& match : matches) {
            outFile << query << '\t' << match.junctionAA << '\t' << match.distance;
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
    trie.SetDeletionCost(config.deletionCost);

    fs::create_directories(config.outputPath);
    std::string outFilePath = config.outputPath + "/results.tsv";

    if (!config.query.empty()) {
        std::vector<AIRREntity> results;
        if (!config.matrixPath.empty()) {
            results = trie.SearchWithMatrix(config.query, config.costRadius);
        } else {
            results = trie.SearchAIRR(config.query, config.maxSubstitution, config.maxInsertion, config.maxDeletion);
        }
        std::unordered_map<std::string, std::vector<AIRREntity>> wrapped{{config.query, results}};
        WriteResults(outFilePath, wrapped);
    }
    else if (!config.inputQueries.empty()) {
        auto queries = LoadQueriesFromFile(config.inputQueries);
        const size_t BATCH_SIZE = 1000;
        bool firstBatch = true;

        for (size_t start = 0; start < queries.size(); start += BATCH_SIZE) {
            size_t end = std::min(queries.size(), start + BATCH_SIZE);
            std::vector<std::string> batchQueries(queries.begin() + start, queries.begin() + end);

            std::unordered_map<std::string, std::vector<AIRREntity>> batchResults;
            if (!config.matrixPath.empty()) {
                batchResults = trie.SearchForAllWithMatrix(batchQueries, config.costRadius);
            } else {
                batchResults = trie.SearchForAll(batchQueries, config.maxSubstitution, config.maxInsertion, config.maxDeletion);
            }

            std::ofstream outFile(outFilePath, firstBatch ? std::ios::out : std::ios::out | std::ios::app);
            if (!outFile.is_open()) {
                std::cerr << "Error: Unable to write to " << outFilePath << std::endl;
                return;
            }

            if (firstBatch) {
                bool hasVGene = false, hasJGene = false;
                for (const auto& [_, matches] : batchResults) {
                    for (const auto& m : matches) {
                        if (!m.vGene.empty()) hasVGene = true;
                        if (!m.jGene.empty()) hasJGene = true;
                        if (hasVGene && hasJGene) break;
                    }
                    if (hasVGene && hasJGene) break;
                }
                outFile << "query\tmatch\tdist";
                if (hasVGene) outFile << "\tv_gene";
                if (hasJGene) outFile << "\tj_gene";
                outFile << '\n';
            }

            for (const auto& [query, matches] : batchResults) {
                for (const auto& m : matches) {
                    outFile << query << '\t'
                            << m.junctionAA << '\t'
                            << m.distance;
                    if (!m.vGene.empty()) outFile << '\t' << m.vGene;
                    if (!m.jGene.empty()) outFile << '\t' << m.jGene;
                    outFile << '\n';
                }
            }

            firstBatch = false;
        }
    }
    else {
        std::cerr << "Error: No query provided.\n";
    }

    std::cout << "SearchAIRR complete. Results saved to: " << outFilePath << std::endl;
}
