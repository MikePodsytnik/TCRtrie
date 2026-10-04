#include "TsvReader.h"

#include <cctype>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
    std::string Trim(std::string s) {
        while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
        std::size_t start = 0;
        while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) ++start;
        return s.substr(start);
    }

    void TrimCR(std::string& s) {
        if (!s.empty() && s.back() == '\r') s.pop_back();
    }

    std::vector<std::string> SplitTSV(const std::string& line) {
        std::vector<std::string> fields;
        std::size_t start = 0;
        while (start <= line.size()) {
            std::size_t end = line.find('\t', start);
            if (end == std::string::npos) end = line.size();
            fields.emplace_back(line.substr(start, end - start));
            start = end + 1;
            if (end == line.size()) break;
        }
        return fields;
    }

    std::string NormalizeSequence(std::string sequence) {
        for (char& c : sequence) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        return sequence;
    }

    bool IsValidSequence(const std::string& sequence) {
        if (sequence.empty()) return false;
        for (char c : sequence) {
            if (c < 'A' || c > 'Z') return false;
        }
        return true;
    }

    int ResolveColumn(const std::unordered_map<std::string, int>& columns,
                      const std::optional<std::string>& explicitName,
                      const std::vector<std::string>& fallbackNames,
                      bool required,
                      const std::string& errorName) {
        if (explicitName) {
            auto it = columns.find(*explicitName);
            if (it == columns.end()) {
                if (required) throw std::runtime_error("Missing required column: " + *explicitName);
                return -1;
            }
            return it->second;
        }
        for (const auto& name : fallbackNames) {
            auto it = columns.find(name);
            if (it != columns.end()) return it->second;
        }
        if (required) throw std::runtime_error("Missing required column: " + errorName);
        return -1;
    }

    int FindColumn(const std::unordered_map<std::string, int>& columns,
                   const std::vector<std::string>& names) {
        for (const auto& name : names) {
            auto it = columns.find(name);
            if (it != columns.end()) return it->second;
        }
        return -1;
    }

    std::string FieldAt(const std::vector<std::string>& fields, int idx) {
        if (idx < 0 || idx >= static_cast<int>(fields.size())) return {};
        return fields[idx];
    }

    std::optional<int> ParseOptionalInt(const std::string& value) {
        std::string trimmed = Trim(value);
        if (trimmed.empty()) return std::nullopt;
        try {
            std::size_t pos = 0;
            int parsed = std::stoi(trimmed, &pos);
            if (pos != trimmed.size()) return std::nullopt;
            return parsed;
        } catch (...) {
            return std::nullopt;
        }
    }

    std::optional<int> ExtractJsonInt(const std::string& json, const std::string& key) {
        const std::string quotedKey = "\"" + key + "\"";
        std::size_t keyPos = json.find(quotedKey);
        if (keyPos == std::string::npos) return std::nullopt;

        std::size_t colon = json.find(':', keyPos + quotedKey.size());
        if (colon == std::string::npos) return std::nullopt;

        std::size_t pos = colon + 1;
        while (pos < json.size() && std::isspace(static_cast<unsigned char>(json[pos]))) ++pos;

        bool quoted = false;
        if (pos < json.size() && json[pos] == '"') {
            quoted = true;
            ++pos;
        }

        std::size_t begin = pos;
        if (pos < json.size() && (json[pos] == '-' || json[pos] == '+')) ++pos;
        while (pos < json.size() && std::isdigit(static_cast<unsigned char>(json[pos]))) ++pos;

        if (begin == pos) return std::nullopt;
        if (quoted && (pos >= json.size() || json[pos] != '"')) return std::nullopt;

        return ParseOptionalInt(json.substr(begin, pos - begin));
    }

    std::pair<std::optional<int>, std::optional<int>> ReadRegionBounds(
        const std::vector<std::string>& fields,
        int vEndCol,
        int jStartCol,
        int cdr3FixCol) {

        std::optional<int> vEnd;
        std::optional<int> jStart;

        if (vEndCol >= 0) vEnd = ParseOptionalInt(FieldAt(fields, vEndCol));
        if (jStartCol >= 0) jStart = ParseOptionalInt(FieldAt(fields, jStartCol));

        if (vEnd && jStart) return {vEnd, jStart};

        if (cdr3FixCol >= 0) {
            const std::string cdr3fix = FieldAt(fields, cdr3FixCol);
            auto jsonVEnd = ExtractJsonInt(cdr3fix, "vEnd");
            auto jsonJStart = ExtractJsonInt(cdr3fix, "jStart");
            if (jsonVEnd && jsonJStart) return {jsonVEnd, jsonJStart};
        }

        return {std::nullopt, std::nullopt};
    }

    std::vector<Region> BuildRegions(std::size_t length, int vEnd, int jStart) {
        if (length == 0) throw std::runtime_error("Cannot build regions for empty CDR3");
        if (vEnd < 1 || jStart <= vEnd || jStart > static_cast<int>(length)) {
            throw std::runtime_error("Invalid 1-based vEnd/jStart markup");
        }
        std::vector<Region> regions;
        regions.reserve(length);
        for (std::size_t i = 0; i < length; ++i) {
            int position = static_cast<int>(i) + 1;
            if (position <= vEnd) regions.push_back(Region::V);
            else if (position < jStart) regions.push_back(Region::NDN);
            else regions.push_back(Region::J);
        }
        return regions;
    }
}

TsvReadResult ReadTsv(const std::string& path, const CliConfig& config) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Failed to open TSV: " + path);

    TsvReadResult result;
    std::string headerLine;
    if (!std::getline(input, headerLine)) throw std::runtime_error("Empty TSV: " + path);

    TrimCR(headerLine);
    result.header = SplitTSV(headerLine);

    std::unordered_map<std::string, int> columns;
    for (int i = 0; i < static_cast<int>(result.header.size()); ++i) columns.emplace(result.header[i], i);

    const int junctionCol = ResolveColumn(columns, config.junctionCol, {"junction_aa", "cdr3"}, true, "junction_aa or cdr3");
    const int vCol = ResolveColumn(columns, config.vCol, {"v_call", "v.segm"}, false, "V gene");
    const int jCol = ResolveColumn(columns, config.jCol, {"j_call", "j.segm"}, false, "J gene");
    const int epitopeCol = ResolveColumn(columns, std::optional<std::string>(config.epitopeCol), {}, false, "epitope");
    const int speciesCol = ResolveColumn(columns, std::optional<std::string>(config.speciesCol), {}, false, "species");
    const int chainCol = ResolveColumn(columns, std::optional<std::string>(config.chainCol), {}, false, "chain");
    const int vEndCol = FindColumn(columns, {config.vEndCol, "vEnd", "v_end"});
    const int jStartCol = FindColumn(columns, {config.jStartCol, "jStart", "j_start"});
    const int cdr3FixCol = FindColumn(columns, {config.cdr3FixCol, "cdr3fix"});

    if (config.matchV && vCol < 0) throw std::runtime_error("Flag --match-v was requested, but V column is missing");
    if (config.matchJ && jCol < 0) throw std::runtime_error("Flag --match-j was requested, but J column is missing");
    if (config.epitope && epitopeCol < 0) throw std::runtime_error("Epitope filter was requested, but epitope column is missing: " + config.epitopeCol);
    if (config.species && speciesCol < 0) throw std::runtime_error("Species filter was requested, but species column is missing: " + config.speciesCol);
    if (config.gene && chainCol < 0) throw std::runtime_error("Gene filter was requested, but chain column is missing: " + config.chainCol);

    std::string line;
    std::size_t dataRowIndex = 0;

    while (std::getline(input, line)) {
        TrimCR(line);
        if (line.empty()) continue;
        auto fields = SplitTSV(line);

        Record record;
        record.rowIndex = dataRowIndex++;
        record.junctionAA = NormalizeSequence(Trim(FieldAt(fields, junctionCol)));
        record.vGene = Trim(FieldAt(fields, vCol));
        record.jGene = Trim(FieldAt(fields, jCol));
        record.epitope = Trim(FieldAt(fields, epitopeCol));
        record.species = Trim(FieldAt(fields, speciesCol));
        record.chain = Trim(FieldAt(fields, chainCol));
        record.rawFields = std::move(fields);

        if (record.junctionAA.empty()) {
            ++result.skippedMissingJunction;
            continue;
        }
        if (!IsValidSequence(record.junctionAA)) {
            ++result.skippedInvalidSequence;
            continue;
        }

        if (config.IsRegionalMatrixMode()) {
            auto [vEnd, jStart] = ReadRegionBounds(record.rawFields, vEndCol, jStartCol, cdr3FixCol);
            if (!vEnd || !jStart) {
                ++result.skippedMissingRegion;
                continue;
            }
            try {
                record.vEnd = *vEnd;
                record.jStart = *jStart;
                record.regions = BuildRegions(record.junctionAA.size(), *vEnd, *jStart);
            } catch (...) {
                ++result.skippedInvalidRegion;
                continue;
            }
        }

        if (config.matchV && record.vGene.empty()) {
            ++result.skippedByFilter;
            continue;
        }
        if (config.matchJ && record.jGene.empty()) {
            ++result.skippedByFilter;
            continue;
        }
        if (config.gene && record.chain != *config.gene) {
            ++result.skippedByFilter;
            continue;
        }
        if (config.species && record.species != *config.species) {
            ++result.skippedByFilter;
            continue;
        }
        if (config.epitope && record.epitope != *config.epitope) {
            ++result.skippedByFilter;
            continue;
        }

        result.records.push_back(std::move(record));
    }

    return result;
}
