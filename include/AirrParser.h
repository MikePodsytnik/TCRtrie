#pragma once

#include <tuple>
#include <string>
#include <string_view>
#include <vector>

struct AIRREntity {
    std::string junctionAA;
    std::string vGene;
    std::string jGene;
    double distance{};

    AIRREntity() = default;
    AIRREntity(std::string_view ja,
               std::string_view v,
               std::string_view j,
               double d)
            : junctionAA(ja)
            , vGene(v)
            , jGene(j)
            , distance(d)
    {}

    auto key() const noexcept { return std::tie(junctionAA, vGene, jGene); }
};

namespace detail {
    template <class T>
    inline void hash_combine(std::size_t& seed, const T& v) {
        seed ^= std::hash<T>{}(v) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    }
}

namespace std {
    template<>
    struct hash<AIRREntity> {
        size_t operator()(const AIRREntity& x) const noexcept {
            size_t seed = 0;
            detail::hash_combine(seed, x.junctionAA);
            detail::hash_combine(seed, x.vGene);
            detail::hash_combine(seed, x.jGene);
            return seed;
        }
    };
}

inline bool operator==(const AIRREntity& a, const AIRREntity& b) noexcept {
    return a.key() == b.key();
}

std::vector<AIRREntity> ParseAIRR(const std::string& filepath);
