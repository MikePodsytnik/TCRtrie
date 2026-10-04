#include "CostMatrix.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
    const std::string kAminoAcids = "ACDEFGHIKLMNPQRSTVWY";

    std::string Trim(const std::string& s) {
        std::size_t left = 0;
        while (left < s.size() && std::isspace(static_cast<unsigned char>(s[left]))) ++left;
        std::size_t right = s.size();
        while (right > left && std::isspace(static_cast<unsigned char>(s[right - 1]))) --right;
        return s.substr(left, right - left);
    }

    void RemoveBom(std::string& s) {
        if (s.size() >= 3 && static_cast<unsigned char>(s[0]) == 0xEF && static_cast<unsigned char>(s[1]) == 0xBB && static_cast<unsigned char>(s[2]) == 0xBF) {
            s.erase(0, 3);
        }
    }

    std::vector<std::string> SplitLine(const std::string& line, const std::string& delimiter) {
        std::vector<std::string> result;
        if (delimiter.empty()) {
            std::istringstream in(line);
            std::string token;
            while (in >> token) result.push_back(token);
            return result;
        }
        std::size_t start = 0;
        while (start <= line.size()) {
            std::size_t pos = line.find(delimiter, start);
            if (pos == std::string::npos) pos = line.size();
            result.push_back(Trim(line.substr(start, pos - start)));
            start = pos + delimiter.size();
            if (pos == line.size()) break;
        }
        return result;
    }

    char ParseLabel(std::string token) {
        RemoveBom(token);
        token = Trim(token);
        if (token.size() != 1) throw std::runtime_error("Invalid matrix label: " + token);
        return static_cast<char>(std::toupper(static_cast<unsigned char>(token[0])));
    }

    void RequireAminoAcidOrGap(char value) {
        if (kAminoAcids.find(value) == std::string::npos && value != '-') {
            throw std::runtime_error(std::string("Unexpected matrix label: ") + value);
        }
    }
}

void CostMatrix::Load(const std::string& path, const std::string& delimiter, float gapFactor) {
    if (gapFactor < 1.0f) throw std::runtime_error("gapFactor must be >= 1.0");

    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open matrix: " + path);

    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("Matrix file is empty: " + path);

    auto header = SplitLine(line, delimiter);
    if (header.empty()) throw std::runtime_error("Matrix header is empty: " + path);
    RemoveBom(header[0]);

    if (!delimiter.empty()) {
        if (!Trim(header[0]).empty()) throw std::runtime_error("Top-left matrix cell must be empty: " + path);
        header.erase(header.begin());
    }

    if (header.size() != 20 && header.size() != 21) {
        throw std::runtime_error("Matrix must contain 20 or 21 columns: " + path);
    }

    std::vector<char> cols;
    for (const auto& token : header) {
        char label = ParseLabel(token);
        RequireAminoAcidOrGap(label);
        cols.push_back(label);
    }

    std::set<char> colSet(cols.begin(), cols.end());
    if (colSet.size() != cols.size()) throw std::runtime_error("Duplicate matrix column labels: " + path);

    for (char aa : kAminoAcids) {
        if (colSet.find(aa) == colSet.end()) throw std::runtime_error("Matrix misses amino acid column: " + std::string(1, aa));
    }

    std::unordered_map<char, std::unordered_map<char, float>> score;
    std::vector<char> rows;

    while (std::getline(input, line)) {
        if (Trim(line).empty()) continue;
        auto parts = SplitLine(line, delimiter);
        if (parts.size() != cols.size() + 1) throw std::runtime_error("Invalid matrix row width: " + path);
        char row = ParseLabel(parts[0]);
        RequireAminoAcidOrGap(row);
        rows.push_back(row);
        for (std::size_t i = 0; i < cols.size(); ++i) {
            std::size_t pos = 0;
            float value = std::stof(Trim(parts[i + 1]), &pos);
            if (pos != Trim(parts[i + 1]).size() || !std::isfinite(value)) throw std::runtime_error("Invalid matrix value: " + path);
            score[row][cols[i]] = value;
        }
    }

    std::set<char> rowSet(rows.begin(), rows.end());
    if (rows.size() != cols.size()) throw std::runtime_error("Matrix is not square: " + path);
    if (rowSet.size() != rows.size()) throw std::runtime_error("Duplicate matrix row labels: " + path);
    if (rowSet != colSet) throw std::runtime_error("Matrix row and column labels do not match: " + path);

    for (char aa : kAminoAcids) {
        float diag = score.at(aa).at(aa);
        for (char other : cols) {
            if (other == aa) continue;
            if (diag <= score.at(aa).at(other) || diag <= score.at(other).at(aa)) {
                throw std::runtime_error("Matrix diagonal must be strictly greater for amino acid: " + std::string(1, aa));
            }
        }
    }

    std::vector<char> labels = cols;
    bool hasGap = colSet.find('-') != colSet.end();

    if (!hasGap) {
        score['-']['-'] = 0.0f;
        for (char aa : kAminoAcids) {
            float minValue = std::numeric_limits<float>::infinity();
            for (char rowAa : kAminoAcids) {
                minValue = std::min(minValue, score.at(rowAa).at(aa));
            }
            float gapScore = minValue < 0.0f ? minValue * gapFactor : minValue / gapFactor;
            score[aa]['-'] = gapScore;
            score['-'][aa] = gapScore;
        }
        labels.push_back('-');
    }

    cost_.clear();
    for (char r : labels) {
        for (char c : labels) {
            float value = 0.0f;
            if (r == '-' && c == '-') {
                value = 0.0f;
            } else if (r == '-') {
                value = score.at(c).at(c) + std::abs(score.at(r).at(c));
            } else if (c == '-') {
                value = score.at(r).at(r) + std::abs(score.at(r).at(c));
            } else {
                value = (score.at(r).at(r) + score.at(c).at(c)) * 0.5f - score.at(r).at(c);
            }
            if (value < 0.0f) throw std::runtime_error("Negative cost after conversion in matrix: " + path);
            cost_[r][c] = value;
        }
    }
}

float CostMatrix::Cost(char row, char column) const {
    char r = static_cast<char>(std::toupper(static_cast<unsigned char>(row)));
    char c = static_cast<char>(std::toupper(static_cast<unsigned char>(column)));
    auto rowIt = cost_.find(r);
    if (rowIt == cost_.end()) throw std::runtime_error(std::string("No matrix row for: ") + r);
    auto colIt = rowIt->second.find(c);
    if (colIt == rowIt->second.end()) throw std::runtime_error(std::string("No matrix column for: ") + c);
    return colIt->second;
}

bool CostMatrix::Loaded() const {
    return !cost_.empty();
}
