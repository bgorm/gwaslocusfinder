#pragma once
#include "core.hpp"
#include <htslib/hts.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
// Earlier pinned releases can crash while cleaning up a failed indexed-reader
// open. 1.25 fixes ownership/cleanup on these error paths.
#if HTS_VERSION < 102500
#error gwaslocusfinder requires HTSlib 1.25 or newer
#endif
#include <zlib.h>
#include <zstd.h>
#include <array>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <unordered_map>

namespace lf {
namespace fs = std::filesystem;
inline bool ends_with(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() && s.compare(s.size()-suffix.size(), suffix.size(), suffix) == 0;
}
inline std::vector<std::string> fields(const std::string& s) {
    std::istringstream in(s);
    std::vector<std::string> out;
    for (std::string v; in >> v;) out.push_back(v);
    return out;
}
// Streaming plain/gzip/BGZF/zstd text, including concatenated zstd frames.
class TextReader {
    std::string path_, pending_;
    size_t cursor_ = 0;
    gzFile gz_ = nullptr;
    FILE* file_ = nullptr;
    ZSTD_DStream* zst_ = nullptr;
    std::array<char, 65536> input_{}, output_{};
    ZSTD_inBuffer in_{nullptr, 0, 0};
    size_t remaining_ = 0;
    bool finished_ = false;
    size_t refill() {
        if (gz_) {
            int n = gzread(gz_, output_.data(), output_.size());
            int err = Z_OK;
            gzerror(gz_, &err);
            require(n >= 0 && (err == Z_OK || err == Z_STREAM_END), "Cannot decompress " + path_);
            return static_cast<size_t>(n);
        }
        while (true) {
            if (in_.pos == in_.size) {
                size_t n = fread(input_.data(), 1, input_.size(), file_);
                require(!ferror(file_), "Cannot read " + path_);
                if (!n) {
                    require(remaining_ == 0, "Truncated zstd file: " + path_);
                    return 0;
                }
                in_ = {input_.data(), n, 0};
            }
            ZSTD_outBuffer out{output_.data(), output_.size(), 0};
            remaining_ = ZSTD_decompressStream(zst_, &out, &in_);
            require(!ZSTD_isError(remaining_), "Invalid zstd file: " + path_);
            if (out.pos) return out.pos;
        }
    }
public:
    explicit TextReader(const std::string& path) : path_(path) {
        if (ends_with(path, ".zst")) {
            file_ = fopen(path.c_str(), "rb");
            require(file_, "Cannot open " + path);
            zst_ = ZSTD_createDStream();
            require(zst_, "Cannot allocate zstd reader");
            remaining_ = ZSTD_initDStream(zst_);
            require(!ZSTD_isError(remaining_), "Cannot initialize zstd reader");
        } else {
            gz_ = gzopen(path.c_str(), "rb");
            require(gz_, "Cannot open " + path);
        }
    }
    TextReader(const TextReader&) = delete;
    ~TextReader() { if (gz_) gzclose(gz_); if (file_) fclose(file_); if (zst_) ZSTD_freeDStream(zst_); }
    bool getline(std::string& line) {
        while (true) {
            auto end = pending_.find('\n', cursor_);
            if (end != std::string::npos) {
                line = pending_.substr(cursor_, end-cursor_); cursor_ = end+1;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
            if (finished_) {
                if (cursor_ == pending_.size()) return false;
                line = pending_.substr(cursor_); cursor_ = pending_.size();
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
            pending_.erase(0, cursor_); cursor_ = 0;
            auto n = refill();
            if (n) pending_.append(output_.data(), n); else finished_ = true;
        }
    }
};
struct RecordDeleter { void operator()(bcf1_t* p) const { bcf_destroy(p); } };
struct HeaderDeleter { void operator()(bcf_hdr_t* p) const { bcf_hdr_destroy(p); } };
using Record = std::shared_ptr<bcf1_t>;
using Header = std::unique_ptr<bcf_hdr_t, HeaderDeleter>;
// Indexed local VCF/BCF reader shared by GWAS and reference inputs.
class IndexedReader {
    bcf_srs_t* sr_ = nullptr;
    int last_rid_ = -1;
    hts_pos_t last_pos_ = -1;
    std::unordered_map<int, bool> visited_;
    std::string path_;
    void reset_order() { last_rid_ = -1; last_pos_ = -1; visited_.clear(); }
public:
    bcf_hdr_t* header = nullptr;
    IndexedReader(const std::string& path, const std::string& region, int threads,
                  bool regions_file = false, int regions_overlap = 0,
                  const std::string& targets = "", bool targets_file = false, int targets_overlap = 0) : path_(path) {
        require(path != "-" && fs::is_regular_file(path),
                "Input must be an indexed, rereadable local VCF/BCF file: " + path);
        sr_ = bcf_sr_init();
        require(sr_, "Cannot allocate HTSlib reader");
        try {
            bcf_sr_set_opt(sr_, BCF_SR_REQUIRE_IDX);
            bcf_sr_set_opt(sr_, BCF_SR_REGIONS_OVERLAP, regions_overlap);
            bcf_sr_set_opt(sr_, BCF_SR_TARGETS_OVERLAP, targets_overlap);
            if (threads > 0) require(bcf_sr_set_threads(sr_, threads) == 0, "Cannot start HTSlib threads");
            if (!region.empty())
                require(bcf_sr_set_regions(sr_, region.c_str(), regions_file) == 0, "Invalid regions: " + region);
            if (!targets.empty())
                require(bcf_sr_set_targets(sr_,targets.c_str(),targets_file,0) == 0, "Invalid targets: " + targets);
            require(bcf_sr_add_reader(sr_, path.c_str()),
                    "Cannot open indexed VCF/BCF " + path + "; BGZF compression and a matching CSI/TBI index are required");
            header = bcf_sr_get_header(sr_, 0);
        } catch (...) { bcf_sr_destroy(sr_); sr_ = nullptr; throw; }
    }
    IndexedReader(const IndexedReader&) = delete;
    virtual ~IndexedReader() { if (sr_) bcf_sr_destroy(sr_); }
    void regions(const std::string& region) {
        require(!region.empty() && bcf_sr_set_regions(sr_,region.c_str(),0) == 0, "Cannot set reference regions");
        reset_order();
    }
    void seek(const std::string& chr, int64_t pos) {
        require(bcf_sr_seek(sr_,chr.c_str(),pos-1) == 0,
                "Cannot seek reference " + path_ + " at " + chr + ":" + std::to_string(pos));
        reset_order();
    }
    bcf1_t* next() {
        if (bcf_sr_next_line(sr_)) {
            auto* r = bcf_sr_get_line(sr_, 0);
            require(!r->errcode, "Malformed VCF/BCF record in " + path_);
            bcf_unpack(r, BCF_UN_STR);
            if (r->rid != last_rid_) {
                require(!visited_.count(r->rid), "Unsorted chromosome blocks in " + path_);
                visited_[r->rid] = true; last_rid_ = r->rid; last_pos_ = -1;
            }
            require(r->pos >= last_pos_, "Unsorted positions in " + path_);
            last_pos_ = r->pos;
            return r;
        }
        require(!sr_->errnum, "VCF/BCF read failed: " + std::string(bcf_sr_strerror(sr_->errnum)));
        return nullptr;
    }
};
class GwasReader : public IndexedReader {
    float* values_ = nullptr;
    int capacity_ = 0;
public:
    GwasReader(const std::string& path, const std::string& region, int threads,
               bool regions_file = false, int regions_overlap = 1,
               const std::string& targets = "", bool targets_file = false, int targets_overlap = 0)
        : IndexedReader(path,region,threads,regions_file,regions_overlap,targets,targets_file,targets_overlap) {
        require(bcf_hdr_nsamples(header) > 0, "GWAS-VCF must have phenotype/sample columns");
        int id = bcf_hdr_id2int(header, BCF_DT_ID, "LP");
        require(id >= 0 && bcf_hdr_idinfo_exists(header, BCF_HL_FMT, id), "Missing FORMAT/LP header");
        require(bcf_hdr_id2type(header, BCF_HL_FMT, id) == BCF_HT_REAL, "FORMAT/LP must have Type=Float");
    }
    ~GwasReader() override { free(values_); }
    void scalars(bcf1_t* r, const char* tag, const std::vector<int>& samples, std::vector<double>& out) {
        out.assign(samples.size(),std::numeric_limits<double>::quiet_NaN());
        int n = bcf_get_format_float(header,r,tag,&values_,&capacity_);
        if (n <= 0) return;
        require(n == bcf_hdr_nsamples(header), std::string("Expected one FORMAT/") + tag + " value per phenotype at " +
                bcf_hdr_id2name(header,r->rid) + ":" + std::to_string(r->pos+1));
        for (size_t i = 0; i < samples.size(); ++i) {
            float v = values_[samples[i]];
            if (!bcf_float_is_missing(v) && !bcf_float_is_vector_end(v)) out[i] = v;
        }
    }
};
} // namespace lf
