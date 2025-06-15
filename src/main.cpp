#include <CLI/CLI.hpp>
#include "TrieInterface.h"
#include <iostream>

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    CLI::App app{"TCRtrie CLI – Approximate TCR sequence search"};

    SearchConfig config;

    app.add_option("-t,--trie", config.inputPath, "Path to AIRR file with sequences")->required();
    app.add_option("-o,--output", config.outputPath, "Path to output folder");

    auto* queryOpt = app.add_option("-q,--query", config.query, "Single query sequence");
    app.add_option("--v-gene", config.vGene, "V-gene to match")->needs(queryOpt);
    app.add_option("--j-gene", config.jGene, "J-gene to match")->needs(queryOpt);
    app.add_option("--input-queries", config.inputQueries, "Path to AIRR file with batch query sequences");

    app.add_option("-s,--sub", config.maxSubstitution, "Allowable number of substitutions");
    app.add_option("-i,--ins", config.maxInsertion, "Allowed number of inserts");
    app.add_option("-d,--del", config.maxDeletion, "Allowed number of deletions");

    auto* matrixOpt = app.add_option("-m,--matrix-search", config.matrixPath, "Path to substitution matrix file");
    app.add_option("-r,--score-radius", config.costRadius, "Score radius for matrix-based search")->needs(matrixOpt);

    app.callback([&]() {
        if (config.query.empty() && config.inputQueries.empty()) {
            throw CLI::ValidationError("No query received");
        }

        if (!config.query.empty() && !config.inputQueries.empty()) {
            throw CLI::ValidationError("Only one of --query or --input-queries must be specified.");
        }

        if (!config.matrixPath.empty() && (config.maxSubstitution >= 0
                                        || config.maxInsertion >= 0
                                        || config.maxDeletion >= 0)) {
            throw CLI::ValidationError("Only one of Levenshtein or Score search must be specified.");
        }

        if (!config.matrixPath.empty() && config.costRadius < 0) {
            throw CLI::ValidationError("--score-radius must be specified with --matrix-search.");
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
        std::cout << "Execution time" << elapsed_time_classic.count() << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "Error during search: " << ex.what() << std::endl;
        return 1;
    }

    return 0;
}