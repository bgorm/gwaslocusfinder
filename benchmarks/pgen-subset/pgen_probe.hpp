#pragma once
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

// Included only in the instrumented builds used for the benchmark.
struct PgenReadProbe {
    using Clock = std::chrono::steady_clock;
    Clock::duration library{}, row_total{};
    uint64_t rows = 0, dosage_rows = 0, dosage_values = 0;
    bool emit_rows = std::getenv("GWASLOCUSFINDER_BENCH_VALUES") != nullptr;

    void observe(Clock::time_point begin, Clock::time_point decoded, Clock::time_point ready,
                 uint32_t row, const std::vector<float>& values, uint32_t dosage_count) {
        library += decoded - begin;
        row_total += ready - begin;
        ++rows;
        dosage_rows += dosage_count != 0;
        dosage_values += dosage_count;
        if (emit_rows &&
            (std::fwrite(&row, sizeof(row), 1, stdout) != 1 ||
             std::fwrite(values.data(), sizeof(float), values.size(), stdout) != values.size()))
            throw std::runtime_error("Cannot write benchmark row stream");
    }
    ~PgenReadProbe() {
        if (!rows) return;
        std::fprintf(stderr, "PGEN_DECODE\t%llu\t%llu\t%llu\t%.9f\t%.9f\n",
                     static_cast<unsigned long long>(rows), static_cast<unsigned long long>(dosage_rows),
                     static_cast<unsigned long long>(dosage_values), std::chrono::duration<double>(library).count(),
                     std::chrono::duration<double>(row_total).count());
    }
};
