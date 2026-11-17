#include "entry.hpp"
#include <htslib/hts.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <strings.h>
#include <stdexcept>
#ifndef LF_BCFTOOLS_VERSION
#error Build with the target bcftools version
#endif
#ifndef LF_HTSLIB_VERSION
#error Build with the target HTSlib version
#endif
extern "C" const char* bcftools_version(void);
// Opaque expression API from bcftools 1.25 filter.h. Symbols resolve from
// the version-checked host, so expressions use bcftools' own evaluator.
extern "C" {
struct _filter_t;
_filter_t* filter_parse(bcf_hdr_t*, const char*);
int filter_status(_filter_t*);
int filter_test(_filter_t*, bcf1_t*, const uint8_t**);
void filter_destroy(_filter_t*);
char* expand_path(char*);
}
namespace {
// The host evaluates the expression. Track its file operands for the output
// collision check using the filename boundaries in bcftools 1.25 filter.c.
std::vector<std::string> filter_input_files(const std::string& expression) {
    auto expanded = [](std::string name) {
        std::unique_ptr<char,decltype(&std::free)> path(expand_path(name.data()),std::free);
        if (!path) throw std::bad_alloc();
        return std::string(path.get());
    };
    std::vector<std::string> paths;
    constexpr const char* spaces = " \t\n\r\f\v";
    auto pos = expression.find_first_not_of(spaces);
    if (pos == std::string::npos) return paths;
    if (!strncasecmp(expression.c_str()+pos,"perl:",5)) {
        const auto end = expression.find(';',pos+5);
        if (end == std::string::npos) return paths;
        paths.push_back(expanded(expression.substr(pos+5,end-pos-5)));
        pos = end+1;
    }
    while ((pos = expression.find_first_of("@\"'",pos)) != std::string::npos) {
        if (expression[pos] != '@') {
            const auto end = expression.find(expression[pos],pos+1);
            if (end == std::string::npos) break;
            pos = end+1;
            continue;
        }
        const bool sample_list = pos && expression[pos-1] == '[';
        if (pos && !sample_list && !std::strchr(" \t\n\r\f\v=!~<>(),&|+-*/%",expression[pos-1])) {
            ++pos;
            continue;
        }
        const auto start = ++pos;
        pos = sample_list ? expression.find(']',start) : expression.find_first_of(" \t\n\r\f\v=!",start);
        auto path = expanded(expression.substr(start,pos-start));
        if (sample_list && !std::filesystem::exists(path)) {
            // FMT/TAG[@file:indices] tries the literal filename first.
            const auto colon = path.rfind(':');
            if (colon != std::string::npos) path.resize(colon);
        }
        paths.push_back(std::move(path));
    }
    return paths;
}
class BcftoolsFilter final : public lf::RecordFilter {
    std::unique_ptr<_filter_t, decltype(&filter_destroy)> filter_{nullptr, filter_destroy};
public:
    BcftoolsFilter(bcf_hdr_t* header, const std::string& expression) {
        filter_.reset(filter_parse(header,expression.c_str()));
        if (!filter_ || filter_status(filter_.get()) != 0)
            throw std::runtime_error("Invalid bcftools filter expression (check syntax and tag names): " + expression);
    }
    bool matches(bcf1_t* record) override { return filter_test(filter_.get(),record,nullptr) != 0; }
};
std::unique_ptr<lf::RecordFilter> make_filter(bcf_hdr_t* header, const std::string& expression) {
    return std::make_unique<BcftoolsFilter>(header,expression);
}
}
#define LF_EXPORT extern "C" __attribute__((visibility("default")))
LF_EXPORT const char* about(void) {
    return "Find significant loci in GWAS-VCF summary statistics using LD clumping and FUMA-style aggregation.\n";
}
LF_EXPORT const char* usage(void) {
    static const std::string text = lf::usage_text(true);
    return text.c_str();
}
LF_EXPORT void version(const char** bcftools, const char** htslib) {
    *bcftools = LF_BCFTOOLS_VERSION; *htslib = LF_HTSLIB_VERSION;
}
LF_EXPORT int run(int argc, char** argv) {
    if (std::strcmp(hts_version(),LF_HTSLIB_VERSION) || std::strcmp(bcftools_version(),LF_BCFTOOLS_VERSION)) {
        std::fprintf(stderr,"gwaslocusfinder: plugin built for bcftools %s / HTSlib %s; rebuild for host %s / %s\n",
                     LF_BCFTOOLS_VERSION,LF_HTSLIB_VERSION,bcftools_version(),hts_version());
        return 1;
    }
    return lf::entry(argc,argv,true,{make_filter,filter_input_files});
}
