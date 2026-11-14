#pragma once
#include <htslib/vcf.h>
#include <memory>
#include <string>
namespace lf {
struct RecordFilter {
    virtual ~RecordFilter() = default;
    virtual bool matches(bcf1_t* record) = 0;
};
using FilterFactory = std::unique_ptr<RecordFilter> (*)(bcf_hdr_t*, const std::string&);
std::string usage_text(bool plugin);
int entry(int argc, char** argv, bool plugin = false, FilterFactory filters = nullptr);
}
