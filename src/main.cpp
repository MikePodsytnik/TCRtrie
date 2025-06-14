#include <CLI/CLI.hpp>
#include "TrieInterface.h"
#include "Trie.h"
#include <random>
#include <iostream>

void load_patterns(const std::string& filename, std::vector<std::string>& patterns) {
    std::ifstream file(filename);
    std::string line;
    while(std::getline(file, line)) {
        patterns.push_back(line);
    }
}

void sample_clonotypes(const std::vector<std::string>& patterns, std::vector<std::string>& clonotypes, size_t sample_size) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, patterns.size() - 1);

    clonotypes.clear();
    for (int i = 0; i < sample_size; ++i) {
        clonotypes.push_back(patterns[dis(gen)]);
    }
    std::cout << "Сэмплировано " << sample_size << " клонотипов.\n";
}

int main(int argc, char** argv) {
    CLI::App app{"TCRtrie CLI – Approximate TCR sequence search"};

    SearchConfig config;

    app.add_option("-i,--input", config.inputPath, "Path to AIRR file with patterns")->required();
    app.add_option("-o,--output", config.outputPath, "Path to output folder");
    app.add_option("--n_edits", config.nEdits, "Allowed number of edits (Levenshtein distance)");

    app.add_option("-q,--query", config.query, "Single query sequence");
    app.add_option("--input_queries", config.inputQueries, "Path to AIRR file with batch query sequences");

    app.add_option("--v_gene", config.vGene, "V-gene to match");
    app.add_option("--j_gene", config.jGene, "J-gene to match");

    app.add_option("--matrix_search", config.matrixPath, "Path to substitution matrix file");
    app.add_option("--score_radius", config.scoreRadius, "Score radius for matrix-based search")->needs("--matrix_search");

    app.callback([&]() {
        if (!config.query.empty() && !config.inputQueries.empty()) {
            throw CLI::ValidationError("Only one of --query or --input_queries should be specified.");
        }

        if (!config.query.empty() && config.matrixPath.empty() && config.nEdits == 0) {
            throw CLI::ValidationError("Specify either --n_edits or --matrix_search for --query.");
        }

        if (!config.matrixPath.empty() && config.scoreRadius == 0) {
            throw CLI::ValidationError("--score_radius must be specified with --matrix_search.");
        }

        if (config.outputPath.empty()) {
            config.outputPath = "./";
        }
    });

    CLI11_PARSE(app, argc, argv);

    try {
        auto start_time = std::chrono::high_resolution_clock::now();
        RunSearch(config);
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed_time_classic = end_time - start_time;
        std::cout << elapsed_time_classic.count() << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "Error during search: " << ex.what() << std::endl;
        return 1;
    }
//    std::vector<std::string> patterns;
//    load_patterns("cdr3.txt", patterns);
//
//    Trie trie("/Users/mihail/CLionProjects/TCRtrie/vdjdb_airr.tsv");
//    std::ofstream classic_results("Classic.csv");
//    std::ofstream heuristic_results("Myers.csv");
//
//    for (int dist = 1; dist < 5; ++dist) {
//        std::vector<std::string> clonotypes;
//        sample_clonotypes(patterns, clonotypes, 200);
//        for (auto clonotype : clonotypes) {
//            auto start_time = std::chrono::high_resolution_clock::now();
//            auto result = trie.Search(clonotype, dist);
//            auto end_time = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double, std::milli> elapsed_time_classic = end_time - start_time;
//            auto start_time_0 = std::chrono::high_resolution_clock::now();
//            auto result_hev = trie.Search(clonotype, dist);
//            auto end_time_0 = std::chrono::high_resolution_clock::now();
//            std::chrono::duration<double, std::milli> elapsed_time_hev = end_time_0 - start_time_0;
//            classic_results << dist << ", " << elapsed_time_classic.count() << '\n';
//            heuristic_results << dist << ", " << elapsed_time_hev.count() << '\n';
//            std::cout << result_hev.size() << ' ' << result.size() << '\n';
//        }
//    }

    return 0;
}