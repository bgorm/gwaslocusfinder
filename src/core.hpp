#pragma once
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace lf {
inline void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}
inline std::string chromosome(std::string s) {
    if (s.rfind("chr", 0) == 0) s.erase(0, 3);
    if (s == "23") return "X";
    if (s == "24") return "Y";
    if (s == "M" || s == "26") return "MT";
    return s;
}
inline bool autosome(const std::string& s) {
    for (int i = 1; i <= 22; ++i) if (s == std::to_string(i)) return true;
    return false;
}
inline bool sequence(const std::string& s) {
    return !s.empty() && s.find_first_not_of("ACGTNacgtn") == std::string::npos;
}
inline std::string key(const std::string& chr, int64_t pos, std::string a, std::string b) {
    for (char& c : a) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    for (char& c : b) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    // Reversing single-base alleles preserves SNP identity for squared LD.
    // Reversing an indel changes insertion to deletion at the same coordinate;
    // these can be distinct events and must not share a lookup key.
    if (a.size() == 1 && b.size() == 1 && b < a) std::swap(a, b);
    return chromosome(chr) + ":" + std::to_string(pos) + ":" + a + ":" + b;
}
// Pairwise-complete unphased dosage correlation; negative values denote missing.
inline double r_squared(const std::vector<float>& a, const std::vector<float>& b, size_t min_n) {
    require(a.size() == b.size(), "LD vectors have different sample counts");
    size_t n = 0;
    double sa = 0, sb = 0;
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i] >= 0 && b[i] >= 0) { ++n; sa += a[i]; sb += b[i]; }
    if (n < min_n) return std::numeric_limits<double>::quiet_NaN();
    const double ma = sa / n, mb = sb / n;
    double va = 0, vb = 0, cov = 0;
    for (size_t i = 0; i < a.size(); ++i) if (a[i] >= 0 && b[i] >= 0) {
        const double x = a[i] - ma, y = b[i] - mb;
        va += x*x; vb += y*y; cov += x*y;
    }
    if (va <= 0 || vb <= 0) return std::numeric_limits<double>::quiet_NaN();
    return std::clamp(cov*cov / (va*vb), 0.0, 1.0);
}
struct DisjointSet {
    std::vector<size_t> parent;
    explicit DisjointSet(size_t n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }
    size_t find(size_t i) {
        while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
        return i;
    }
    void join(size_t a, size_t b) { a = find(a); b = find(b); if (a != b) parent[b] = a; }
};
struct Interval { std::string chr; int64_t start, end; size_t item; };
inline void merge_intervals(std::vector<Interval> intervals, int64_t gap, DisjointSet& dsu) {
    std::sort(intervals.begin(), intervals.end(), [](const auto& a, const auto& b) {
        return std::tie(a.chr, a.start, a.end, a.item) < std::tie(b.chr, b.start, b.end, b.item);
    });
    if (intervals.empty()) return;
    auto active = intervals.front();
    for (size_t i = 1; i < intervals.size(); ++i) {
        const auto& next = intervals[i];
        if (next.chr == active.chr && next.start - active.end <= gap) {
            dsu.join(active.item, next.item);
            active.end = std::max(active.end, next.end);
        } else active = next;
    }
}
} // namespace lf
