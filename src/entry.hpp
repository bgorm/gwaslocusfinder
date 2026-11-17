#pragma once
#include <htslib/vcf.h>
#include <memory>
#include <string>
#include <vector>
namespace lf {
struct RecordFilter {
    virtual ~RecordFilter() = default;
    virtual bool matches(bcf1_t* record) = 0;
};
using FilterFactory = std::unique_ptr<RecordFilter> (*)(bcf_hdr_t*, const std::string&);
struct FilterApi {
    FilterFactory compile = nullptr;
    std::vector<std::string> (*input_files)(const std::string&) = nullptr;
};
std::string usage_text(bool plugin);
int entry(int argc, char** argv, bool plugin = false, FilterApi filters = {});
}
