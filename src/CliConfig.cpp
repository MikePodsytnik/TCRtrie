#include "CliConfig.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
    int ParseInt(const std::string& value, const std::string& flag) {
        try {
            std::size_t pos = 0;
            int parsed = std::stoi(value, &pos);
            if (pos != value.size()) throw std::runtime_error("");
            return parsed;
        } catch (...) {
            throw std::runtime_error("Invalid integer for " + flag + ": " + value);
        }
    }

    float ParseFloat(const std::string& value, const std::string& flag) {
        try {
            std::size_t pos = 0;
            float parsed = std::stof(value, &pos);
            if (pos != value.size()) throw std::runtime_error("");
            return parsed;
        } catch (...) {
            throw std::runtime_error("Invalid float for " + flag + ": " + value);
        }
    }

    bool AnyRegionalMatrix(const CliConfig& config) {
        return config.matrixVPath || config.matrixNDNPath || config.matrixJPath;
    }
}

bool CliConfig::IsMatrixMode() const {
    return matrixPath.has_value();
}

bool CliConfig::IsRegionalMatrixMode() const {
    return matrixVPath.has_value() && matrixNDNPath.has_value() && matrixJPath.has_value();
}

void PrintUsage() {
    std::cerr
        << "Usage: tcrtrie <query.tsv> <target.tsv> [options]\n"
        << "Options:\n"
        << "  --out <path>\n"
        << "  --max-sub <int>\n"
        << "  --max-ins <int>\n"
        << "  --max-del <int>\n"
        << "  --max-edits <int>\n"
        << "  --matrix-path <path>\n"
        << "  --matrix-v <path>\n"
        << "  --matrix-ndn <path>\n"
        << "  --matrix-n <path>\n"
        << "  --matrix-j <path>\n"
        << "  --max-cost <float>\n"
        << "  --match-v\n"
        << "  --match-j\n"
        << "  --align\n"
        << "  --gene <value>\n"
        << "  --species <value>\n"
        << "  --epitope <value>\n"
        << "  --threads <int>\n"
        << "  --junction-col <name>\n"
        << "  --v-col <name>\n"
        << "  --j-col <name>\n"
        << "  --epitope-col <name>\n"
        << "  --species-col <name>\n"
        << "  --chain-col <name>\n"
        << "  --v-end-col <name>\n"
        << "  --j-start-col <name>\n"
        << "  --cdr3fix-col <name>\n"
        << "  --help, -h\n";
}

CliConfig ParseCli(int argc, char** argv) {
    if (argc > 1) {
        const std::string firstArg = argv[1];
        if (firstArg == "--help" || firstArg == "-h") {
            PrintUsage();
            std::exit(0);
        }
    }

    if (argc < 3) {
        PrintUsage();
        throw std::runtime_error("Two positional TSV inputs are required");
    }

    CliConfig config;
    config.queryPath = argv[1];
    config.targetPath = argv[2];

    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];

        auto requireValue = [&](const std::string& flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + flag);
            return argv[++i];
        };

        if (arg == "--out") config.outPath = requireValue(arg);
        else if (arg == "--max-sub") config.maxSub = ParseInt(requireValue(arg), arg);
        else if (arg == "--max-ins") config.maxIns = ParseInt(requireValue(arg), arg);
        else if (arg == "--max-del") config.maxDel = ParseInt(requireValue(arg), arg);
        else if (arg == "--max-edits") config.maxEdits = ParseInt(requireValue(arg), arg);
        else if (arg == "--matrix-path") config.matrixPath = requireValue(arg);
        else if (arg == "--matrix-v") config.matrixVPath = requireValue(arg);
        else if (arg == "--matrix-ndn" || arg == "--matrix-n" || arg == "--matrix-d") config.matrixNDNPath = requireValue(arg);
        else if (arg == "--matrix-j") config.matrixJPath = requireValue(arg);
        else if (arg == "--max-cost") config.maxCost = ParseFloat(requireValue(arg), arg);
        else if (arg == "--match-v") config.matchV = true;
        else if (arg == "--match-j") config.matchJ = true;
        else if (arg == "--align") config.writeAlignment = true;
        else if (arg == "--gene") config.gene = requireValue(arg);
        else if (arg == "--species") config.species = requireValue(arg);
        else if (arg == "--epitope") config.epitope = requireValue(arg);
        else if (arg == "--threads") config.threads = ParseInt(requireValue(arg), arg);
        else if (arg == "--junction-col") config.junctionCol = requireValue(arg);
        else if (arg == "--v-col") config.vCol = requireValue(arg);
        else if (arg == "--j-col") config.jCol = requireValue(arg);
        else if (arg == "--epitope-col") config.epitopeCol = requireValue(arg);
        else if (arg == "--species-col") config.speciesCol = requireValue(arg);
        else if (arg == "--chain-col") config.chainCol = requireValue(arg);
        else if (arg == "--v-end-col") config.vEndCol = requireValue(arg);
        else if (arg == "--j-start-col") config.jStartCol = requireValue(arg);
        else if (arg == "--cdr3fix-col") config.cdr3FixCol = requireValue(arg);
        else if (arg == "--help" || arg == "-h") {
            PrintUsage();
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + arg);
        }
    }

    if (config.threads <= 0) throw std::runtime_error("--threads must be positive");
    if (config.maxSub < 0 || config.maxIns < 0 || config.maxDel < 0) throw std::runtime_error("Edit limits must be non-negative");
    if (config.maxEdits < -1) throw std::runtime_error("--max-edits must be non-negative");
    if (config.maxCost < 0.0f) throw std::runtime_error("--max-cost must be non-negative");
    if (config.maxEdits == -1) config.maxEdits = config.maxSub + config.maxIns + config.maxDel;

    if (config.matrixPath && AnyRegionalMatrix(config)) {
        throw std::runtime_error("Use either --matrix-path or regional matrices, not both");
    }

    if (AnyRegionalMatrix(config) && !config.IsRegionalMatrixMode()) {
        throw std::runtime_error("Regional matrix mode requires --matrix-v, --matrix-ndn and --matrix-j");
    }

    return config;
}
