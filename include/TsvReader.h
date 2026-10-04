#pragma once

#include "CliConfig.h"
#include "Record.h"

#include <string>
#include <vector>

struct TsvReadResult {
    std::vector<std::string> header;
    std::vector<Record> records;
    std::size_t skippedMissingJunction = 0;
    std::size_t skippedInvalidSequence = 0;
    std::size_t skippedMissingRegion = 0;
    std::size_t skippedInvalidRegion = 0;
    std::size_t skippedByFilter = 0;
};

TsvReadResult ReadTsv(const std::string& path, const CliConfig& config);
