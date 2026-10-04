#pragma once

#include <optional>
#include <string>

struct CliConfig {
    std::string queryPath;
    std::string targetPath;
    std::string outPath = "match_result.tsv";

    int maxSub = 1;
    int maxIns = 0;
    int maxDel = 0;
    int maxEdits = -1;

    std::optional<std::string> matrixPath = std::nullopt;
    std::optional<std::string> matrixVPath = std::nullopt;
    std::optional<std::string> matrixNDNPath = std::nullopt;
    std::optional<std::string> matrixJPath = std::nullopt;
    float maxCost = 6.0f;

    bool matchV = false;
    bool matchJ = false;
    bool writeAlignment = false;

    std::optional<std::string> gene = std::nullopt;
    std::optional<std::string> species = std::nullopt;
    std::optional<std::string> epitope = std::nullopt;

    int threads = 4;

    std::optional<std::string> junctionCol = std::nullopt;
    std::optional<std::string> vCol = std::nullopt;
    std::optional<std::string> jCol = std::nullopt;
    std::string epitopeCol = "antigen.epitope";
    std::string speciesCol = "species";
    std::string chainCol = "gene";
    std::string vEndCol = "v.end";
    std::string jStartCol = "j.start";
    std::string cdr3FixCol = "cdr3fix";

    bool IsMatrixMode() const;
    bool IsRegionalMatrixMode() const;
};

CliConfig ParseCli(int argc, char** argv);
void PrintUsage();
