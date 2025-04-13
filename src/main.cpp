#include <CLI/CLI.hpp>
#include "TrieInterface.h"
#include <iostream>

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
        RunSearch(config);
    } catch (const std::exception& ex) {
        std::cerr << "Error during search: " << ex.what() << std::endl;
        return 1;
    }

    return 0;
}