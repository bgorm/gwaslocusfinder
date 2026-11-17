#pragma once
#include "storage.hpp"
#ifndef LF_WITH_PGEN
#define LF_WITH_PGEN 1
#endif
#if LF_WITH_PGEN
#include <pgenlib_read.h>
#endif
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <unordered_set>

namespace lf {
#if LF_WITH_PGEN
struct AlignedBuffer {
    void* ptr = nullptr;
    AlignedBuffer() = default;
    AlignedBuffer(const AlignedBuffer&) = delete;
    ~AlignedBuffer() { if (ptr) plink2::aligned_free(ptr); }
    void allocate(size_t bytes) {
        require(!ptr, "Internal error: buffer allocated twice");
        if (bytes) require(!plink2::cachealigned_malloc(bytes, &ptr), "Out of memory allocating PGEN buffers");
    }
    template<class T> T* as() { return static_cast<T*>(ptr); }
};
#endif
using Windows = std::map<std::string, std::vector<int64_t>>;
inline bool near_signal(const Windows& windows, const std::string& chr, int64_t pos, int64_t radius) {
    auto found = windows.find(chr);
    if (found == windows.end()) return false;
    const auto& v = found->second;
    auto it = std::lower_bound(v.begin(), v.end(), pos-radius);
    return it != v.end() && *it <= pos+radius;
}
inline std::vector<std::pair<std::string, std::string>> read_samples(const std::string& path) {
    TextReader in(path);
    std::string line;
    std::vector<std::pair<std::string, std::string>> out;
    std::set<std::pair<std::string, std::string>> seen;
    int iid = -1, fid = -1;
    size_t cols = 0;
    while (in.getline(line)) {
        if (line.empty() || line.rfind("##", 0) == 0) continue;
        auto f = fields(line);
        if (f.empty()) continue;
        if (line[0] == '#') {
            if (iid >= 0) continue;
            f[0].erase(0, 1);
            for (size_t i = 0; i < f.size(); ++i) {
                if (f[i] == "IID") iid = i;
                if (f[i] == "FID") fid = i;
            }
            cols = f.size();
            continue;
        }
        require(iid >= 0 && f.size() == cols, "Expected PSAM header with IID and consistent columns: " + path);
        auto id = std::make_pair(fid < 0 ? "0" : f[fid], f[iid]);
        require(seen.insert(id).second, "Duplicate FID/IID in " + path);
        out.push_back(id);
    }
    require(!out.empty(), "No samples in " + path);
    return out;
}
inline std::vector<uint32_t> select_samples(const std::vector<std::pair<std::string, std::string>>& all,
                                           const std::string& keep) {
    std::vector<uint32_t> selected;
    if (keep.empty()) {
        for (uint32_t i = 0; i < all.size(); ++i) selected.push_back(i);
        return selected;
    }
    TextReader in(keep);
    std::string line;
    std::set<std::pair<std::string, std::string>> pairs;
    std::set<std::string> iids;
    while (in.getline(line)) {
        auto f = fields(line);
        if (f.empty() || f[0][0] == '#') continue;
        require(f.size() == 1 || f.size() == 2, "--keep expects IID or FID IID per line");
        if (f.size() == 1) iids.insert(f[0]); else pairs.emplace(f[0], f[1]);
    }
    require(!pairs.empty() || !iids.empty(), "Empty --keep file");
    std::map<std::string, size_t> counts;
    for (const auto& id : all) ++counts[id.second];
    for (const auto& iid : iids) require(counts[iid] == 1, "--keep IID missing or ambiguous: " + iid);
    size_t found_pairs = 0;
    for (uint32_t i = 0; i < all.size(); ++i) {
        bool match_pair = pairs.count(all[i]);
        found_pairs += match_pair;
        if (match_pair || iids.count(all[i].second)) selected.push_back(i);
    }
    require(found_pairs == pairs.size(), "Some --keep FID/IID pairs are absent from the reference samples");
    return selected;
}
struct GenotypeSource {
    std::string prefix;
    uint32_t samples = 0;
    virtual ~GenotypeSource() = default;
    virtual std::vector<float> read(const RefVariant& variant, const std::vector<uint32_t>& subset, bool hardcalls) = 0;
};
#if LF_WITH_PGEN
struct PgenSource : GenotypeSource {
    uint32_t raw_count = 0, max_alleles = 2;
    std::vector<uintptr_t> offsets;
    plink2::PgenFileInfo info{};
    plink2::PgenReader reader{};
    AlignedBuffer info_mem, reader_mem, nonref, geno, dosage_present, dosage_main;
    bool opened = false;
    PgenSource() { plink2::PreinitPgfi(&info); plink2::PreinitPgr(&reader); }
    PgenSource(const PgenSource&) = delete;
    ~PgenSource() override {
        plink2::PglErr err = plink2::kPglRetSuccess;
        plink2::CleanupPgr(&reader, &err);
        plink2::CleanupPgfi(&info, &err);
    }
    void open() {
        if (opened) return;
        using namespace plink2;
        char error[kPglErrstrBufBlen] = {};
        PgenHeaderCtrl ctrl;
        uintptr_t n_info = 0, n_reader = 0;
        uint32_t max_record = 0;
        const auto path = prefix + ".pgen";
        require(PgfiInitPhase1(path.c_str(), nullptr, raw_count, samples, &ctrl, &info, &n_info, error)
                == kPglRetSuccess, "PGEN initialization failed: " + path + ": " + error);
        if (!offsets.empty()) {
            info.allele_idx_offsets = offsets.data();
            info.max_allele_ct = max_alleles;
        }
        if ((ctrl & 0x30) && offsets.empty()) {
            offsets.resize(static_cast<size_t>(raw_count)+1);
            for (size_t i = 0; i < offsets.size(); ++i) offsets[i] = 2*i;
            info.allele_idx_offsets = offsets.data();
        }
        if ((ctrl & 0xc0) == 0xc0) {
            nonref.allocate(((static_cast<size_t>(raw_count)+511)/512)*64);
            info.nonref_flags = nonref.as<uintptr_t>();
        }
        info_mem.allocate(n_info*64);
        require(PgfiInitPhase2(ctrl, 1, 0, 0, 0, raw_count, &max_record, &info,
                              info_mem.as<unsigned char>(), &n_reader, error) == kPglRetSuccess,
                "PGEN metadata initialization failed: " + path + ": " + error);
        require(!(info.gflags & kfPgenGlobalMultiallelicHardcallFound) || !offsets.empty(),
                "PGEN/PVAR mismatch: multiallelic PGEN has biallelic PVAR: " + path);
        reader_mem.allocate(n_reader*64);
        require(PgrInit(path.c_str(), max_record, &info, &reader, reader_mem.as<unsigned char>())
                == kPglRetSuccess, "Cannot initialize PGEN reader: " + path);
        geno.allocate(((static_cast<size_t>(samples)+255)/256)*64);
        dosage_present.allocate(((static_cast<size_t>(samples)+511)/512)*64);
        dosage_main.allocate(((static_cast<size_t>(samples)*2+63)/64)*64);
        opened = true;
    }
    std::vector<float> read(const RefVariant& variant, const std::vector<uint32_t>& subset, bool hardcalls) override {
        open();
        const uint32_t row = variant.row;
        using namespace plink2;
        PgrSampleSubsetIndex pssi;
        PgrClearSampleSubsetIndex(nullptr, &pssi);
        uint32_t dosage_ct = 0;
        auto* g = geno.as<uintptr_t>();
        auto* present = dosage_present.as<uintptr_t>();
        auto* values = dosage_main.as<uint16_t>();
        PglErr err = hardcalls ? PgrGet(nullptr, pssi, samples, row, &reader, g)
            : PgrGetD(nullptr, pssi, samples, row, &reader, g, present, values, &dosage_ct);
        require(err == kPglRetSuccess, "PGEN read failed: " + prefix + " variant row " + std::to_string(row+1));
        std::vector<float> out(subset.size(), -9);
        size_t dst = 0, d = 0;
        constexpr size_t bits = sizeof(uintptr_t)*8;
        for (uint32_t i = 0; i < samples; ++i) {
            unsigned code = (g[i/(bits/2)] >> (2*(i%(bits/2)))) & 3;
            float value = code == 3 ? -9.0f : static_cast<float>(code);
            if (dosage_ct && ((present[i/bits] >> (i%bits)) & 1)) value = values[d++] / 16384.0f;
            if (dst < subset.size() && subset[dst] == i) out[dst++] = value;
        }
        require(d == dosage_ct, "Inconsistent PGEN dosage count");
        return out;
    }
};
#endif
// Reference VCF records are located through their index on cache misses.
// Metadata contains only variants in the signal windows, never full GT arrays.
struct VcfSource : GenotypeSource {
    int cursor_rid = -1;
    int64_t cursor_pos = -1, cursor_window_end = -1;
    std::map<int,std::vector<std::pair<int64_t,int64_t>>> windows;
    std::unique_ptr<IndexedReader> reader;
    float* ds = nullptr;
    int ds_capacity = 0;
    float* hds = nullptr;
    int hds_capacity = 0;
    int32_t* gt = nullptr;
    int gt_capacity = 0;
    explicit VcfSource(const std::string& path, int threads) {
        prefix = path;
        reader = std::make_unique<IndexedReader>(path,"",threads);
        samples = bcf_hdr_nsamples(reader->header);
        require(samples > 0, "Reference VCF/BCF has no samples: " + path);
    }
    ~VcfSource() override { free(ds); free(hds); free(gt); }
    std::vector<float> read(const RefVariant& variant, const std::vector<uint32_t>& subset, bool hardcalls) override {
        const auto& loc = variant;
        const char* raw_chr = bcf_hdr_id2name(reader->header,loc.row);
        require(raw_chr,"Invalid reference VCF contig ID");
        // Nearby, strictly forward requests can share decompression work.
        // Equal-position requests seek again: key order may differ from file
        // order, and the first duplicate must always win.
        const bool forward = cursor_rid == static_cast<int>(loc.row) &&
            loc.pos > cursor_pos && loc.pos <= cursor_window_end && loc.pos-cursor_pos <= 65536;
        auto find = [&](size_t budget) -> bcf1_t* {
            while (budget--) {
                auto* r = reader->next();
                if (!r || r->rid != static_cast<int>(loc.row) || r->pos+1 > loc.pos) break;
                if (r->pos+1 == loc.pos && r->n_allele == 2 &&
                    key(raw_chr,loc.pos,r->d.allele[0],r->d.allele[1]) == loc.key) return r;
            }
            return nullptr;
        };
        bcf1_t* rec = forward ? find(256) : nullptr;
        if (!rec) { reader->seek(raw_chr,loc.pos); rec = find(SIZE_MAX); }
        require(rec, "Reference record unavailable after indexed seek (check index): " + loc.key);
        cursor_rid = loc.row; cursor_pos = loc.pos;
        const auto& intervals = windows.at(loc.row);
        auto interval = std::upper_bound(intervals.begin(),intervals.end(),loc.pos,
                                         [](int64_t pos,const auto& span) { return pos < span.first; });
        require(interval != intervals.begin() && loc.pos <= std::prev(interval)->second,"Reference outside signal windows");
        cursor_window_end = std::prev(interval)->second;
        int nds = -1, nhds = -1, ngt = bcf_get_genotypes(reader->header,rec,&gt,&gt_capacity);
        if (!hardcalls) {
            int id = bcf_hdr_id2int(reader->header,BCF_DT_ID,"DS");
            if (id >= 0 && bcf_hdr_idinfo_exists(reader->header,BCF_HL_FMT,id))
                require(bcf_hdr_id2type(reader->header,BCF_HL_FMT,id) == BCF_HT_REAL,
                        "Reference FORMAT/DS must have Type=Float");
            nds = bcf_get_format_float(reader->header,rec,"DS",&ds,&ds_capacity);
            id = bcf_hdr_id2int(reader->header,BCF_DT_ID,"HDS");
            if (id >= 0 && bcf_hdr_idinfo_exists(reader->header,BCF_HL_FMT,id))
                require(bcf_hdr_id2type(reader->header,BCF_HL_FMT,id) == BCF_HT_REAL,
                        "Reference FORMAT/HDS must have Type=Float");
            nhds = bcf_get_format_float(reader->header,rec,"HDS",&hds,&hds_capacity);
            require(nhds <= 0 || (nhds % samples == 0 && nhds <= 2*static_cast<int>(samples)),
                    "Reference HDS requires two haplotype dosages per sample: " + loc.key);
            require(nds <= 0 || nds == static_cast<int>(samples), "Reference DS must be scalar per sample: " + loc.key);
        }
        require(ngt > 0 || nds > 0 || nhds > 0, "Reference record requires " + std::string(hardcalls ? "GT" : "DS, HDS or GT") + ": " + loc.key);
        require(ngt <= 0 || ngt % samples == 0, "Invalid reference GT cardinality: " + loc.key);
        int stride = ngt > 0 ? ngt/samples : 0;
        std::vector<float> values(subset.size(),-9);
        for (size_t i = 0; i < subset.size(); ++i) {
            uint32_t s = subset[i];
            int ploidy = 0, called = 0, alt_count = 0;
            for (int j = 0; j < stride; ++j) {
                int32_t g = gt[s*stride+j];
                if (g == bcf_int32_vector_end) break;
                ++ploidy;
                if (bcf_gt_is_missing(g)) continue;
                int allele = bcf_gt_allele(g);
                require(allele == 0 || allele == 1, "Invalid biallelic reference GT: " + loc.key);
                alt_count += allele; ++called;
            }
            require(!called || ploidy == 2,
                    "Reference requires diploid GT, including on chromosome X: " + loc.key +
                    "; use bcftools +fixploidy -- -f 2 to convert GT (DS/HDS are unchanged; "
                    "use --hardcalls if their encoding is incompatible)");
            if (called == 2 && ploidy == 2) values[i] = alt_count;
            if (nhds > 0) {
                bool complete = true;
                double dosage = 0;
                const int hstride = nhds/samples;
                int count = 0;
                bool has_value = false;
                for (int h = 0; h < hstride; ++h) {
                    float d = hds[s*hstride+h];
                    if (bcf_float_is_vector_end(d)) break;
                    ++count;
                    if (bcf_float_is_missing(d)) { complete = false; continue; }
                    require(std::isfinite(d) && d >= 0 && d <= 1, "Reference HDS must be in [0,1]: " + loc.key);
                    has_value = true; dosage += d;
                }
                // A lone '.' is legal missing data, not a haploid dosage.
                require(!has_value || count == 2,
                        "Reference HDS requires two haplotype dosages per sample: " + loc.key);
                if (complete && count == 2) values[i] = dosage;
            }
            if (nds > 0 && !bcf_float_is_missing(ds[s]) && !bcf_float_is_vector_end(ds[s])) {
                require(std::isfinite(ds[s]) && ds[s] >= 0 && ds[s] <= 2,
                        "Reference DS must be in [0,2]: " + loc.key);
                values[i] = ds[s];
            }
        }
        return values;
    }
};
inline std::string signal_regions(bcf_hdr_t* hdr, const Windows& windows, int64_t radius, ReferenceStore* loaded = nullptr) {
    std::ostringstream out;
    bool first = true;
    for (int rid = 0; rid < hdr->n[BCF_DT_CTG]; ++rid) {
        const char* name = bcf_hdr_id2name(hdr,rid);
        if (!name) continue;
        const auto chr = chromosome(name);
        auto found = windows.find(chr);
        if (found == windows.end()) continue;
        int64_t start = 0, end = -1;
        auto append = [&]() {
            if (end < 0 || (loaded && !loaded->overlaps(chr,start,end))) return;
            if (!first) out << ',';
            first = false; out << name << ':' << start << '-' << end;
        };
        for (int64_t pos : found->second) {
            int64_t left = std::max<int64_t>(1,pos-radius), right = pos+radius;
            if (end >= 0 && left <= end+1) end = std::max(end,right);
            else { append(); start = left; end = right; }
        }
        append();
    }
    return out.str();
}
struct ReferenceFilters {
    double maf = 0, call_rate = 0;
    bool passes(const RefQC& qc) const {
        return qc.maf > 0 && qc.maf >= maf && qc.call_rate >= call_rate && qc.variance > 0;
    }
};
struct ReferenceQCCounts {
    size_t evaluated = 0, passed = 0;
    size_t filtered() const { return evaluated-passed; }
};
class Reference {
    std::list<size_t> lru_;
    struct CacheEntry {
        std::shared_ptr<std::vector<float>> data;
        std::list<size_t>::iterator position;
    };
    std::unordered_map<size_t, CacheEntry> cache_;
    size_t cache_limit_;
    bool hardcalls_;
    ReferenceFilters filters_;
public:
    std::vector<std::unique_ptr<GenotypeSource>> shards;
    ScratchDirectory scratch;
    ReferenceStore store;
    std::vector<uint32_t> subset;
    size_t multiallelic_skipped = 0, duplicates_skipped = 0, genotype_reads = 0;
    ReferenceQCCounts qc_counts;
    Reference(const std::vector<std::string>& prefixes, const Windows& windows, int64_t radius,
              const std::string& keep, size_t cache_mb, bool hardcalls, bool vcf_input = false, int threads = 1,
              const std::string& temp_dir = "", size_t phenotypes = 1, ReferenceFilters filters = {})
        : cache_limit_(0), hardcalls_(hardcalls), filters_(filters), scratch(temp_dir), store(scratch,phenotypes) {
        std::vector<std::pair<std::string, std::string>> samples;
        std::set<std::string> unique_prefixes;
        for (auto prefix : prefixes) {
            if (vcf_input) {
                prefix = fs::canonical(prefix).string();
                require(unique_prefixes.insert(prefix).second, "Repeated VCF reference shard: " + prefix);
                auto shard = std::make_unique<VcfSource>(prefix,threads);
                auto* hdr = shard->reader->header;
                std::vector<std::pair<std::string,std::string>> current_samples;
                for (int i = 0; i < bcf_hdr_nsamples(hdr); ++i) current_samples.emplace_back("0",hdr->samples[i]);
                if (shards.empty()) {
                    samples = std::move(current_samples); subset = select_samples(samples,keep);
                } else require(samples == current_samples, "VCF reference shards must have identical sample order: " + prefix);
                load_vcf_variants(*shard,windows,radius,threads);
                shards.push_back(std::move(shard));
                continue;
            }

#if LF_WITH_PGEN
            if (ends_with(prefix, ".pgen")) prefix.resize(prefix.size()-5);
            prefix = fs::weakly_canonical(prefix).string();
            require(unique_prefixes.insert(prefix).second, "Repeated PGEN shard: " + prefix);
            require(fs::exists(prefix+".pgen"), "Missing " + prefix + ".pgen");
            auto shard = std::make_unique<PgenSource>();
            shard->prefix = prefix;
            auto current_samples = read_samples(prefix+".psam");
            if (shards.empty()) {
                samples = std::move(current_samples);
                subset = select_samples(samples, keep);
            } else require(samples == current_samples, "PGEN shards must have identical FID/IID sample order: " + prefix);
            require(samples.size() < UINT32_MAX, "Too many reference samples");
            shard->samples = samples.size();
            load_pgen_variants(*shard,windows,radius);
            shards.push_back(std::move(shard));
#else
            throw std::runtime_error("PGEN support is disabled; rebuild with WITH_PGEN=1 or use --vcf-ref");
#endif
        }
        require(subset.size() >= 3, "Reference sample subset needs at least 3 samples");
        cache_limit_ = std::max<size_t>(1, cache_mb*1024*1024/(subset.size()*sizeof(float)));
        store.build([&](const RefVariant& first, const RefVariant& duplicate) {
            ++duplicates_skipped;
            if (duplicates_skipped <= 10) {
                std::cerr << "Warning: Duplicate/overlapping reference variant: " << duplicate.key
                          << "; keeping first occurrence (" << shards[first.shard]->prefix;
                if (!vcf_input) std::cerr << ", variant row " << first.row+1;
                std::cerr << "), skipping " << shards[duplicate.shard]->prefix;
                if (!vcf_input) std::cerr << ", variant row " << duplicate.row+1;
                std::cerr << '\n';
            }
        });
    }
    std::shared_ptr<std::vector<float>> genotype(size_t i) {
        auto hit = cache_.find(i);
        if (hit != cache_.end()) {
            lru_.splice(lru_.begin(), lru_, hit->second.position);
            return hit->second.data;
        }
        const auto v = store.get(i);
        auto data = std::make_shared<std::vector<float>>(shards[v.shard]->read(v, subset, hardcalls_));
        ++genotype_reads;
        if (!store.qc(i)) {
            double sum = 0, ss = 0;
            size_t n = 0;
            for (float x : *data) if (x >= 0) { sum += x; ++n; }
            double mean = n ? sum/n : 0;
            for (float x : *data) if (x >= 0) ss += (x-mean)*(x-mean);
            const RefQC qc{n ? std::min(mean/2, 1-mean/2) : 0, static_cast<double>(n)/data->size(), ss};
            store.set_qc(i,qc);
            // Persisted QC marks each variant once across cache evictions,
            // repeated LD windows, overlapping shards, and phenotypes.
            ++qc_counts.evaluated;
            qc_counts.passed += filters_.passes(qc);
        }
        while (cache_.size() >= cache_limit_) { cache_.erase(lru_.back()); lru_.pop_back(); }
        lru_.push_front(i);
        cache_[i] = {data, lru_.begin()};
        return data;
    }
    size_t qc_not_evaluated() const { return store.size()-qc_counts.evaluated; }
    bool passes(size_t i) {
        auto qc = store.qc(i);
        if (!qc) { genotype(i); qc = store.qc(i); }
        require(qc.has_value(),"Missing reference QC");
        return filters_.passes(*qc);
    }
private:
    void add_variant(RefVariant variant) {
        variant.shard = shards.size();
        store.add(std::move(variant));
    }
    void load_vcf_variants(VcfSource& shard, const Windows& windows, int64_t radius, int threads) {
        const auto& prefix = shard.prefix;
        const auto regions = signal_regions(shard.reader->header,windows,radius);
        if (regions.empty()) return;
        for (int rid = 0; rid < shard.reader->header->n[BCF_DT_CTG]; ++rid) {
            auto found = windows.find(chromosome(bcf_hdr_id2name(shard.reader->header,rid)));
            if (found == windows.end()) continue;
            auto& intervals = shard.windows[rid];
            for (auto pos : found->second) {
                const int64_t start = std::max<int64_t>(1,pos-radius), end = pos+radius;
                if (!intervals.empty() && start <= intervals.back().second+1) intervals.back().second = end;
                else intervals.emplace_back(start,end);
            }
        }
        IndexedReader scan(prefix,regions,threads,false,0,"",false,0,BCF_UN_INFO);
        while (auto* r = scan.next()) {
            if (r->n_allele != 2) { ++multiallelic_skipped; continue; }
            if (!sequence(r->d.allele[0]) || !sequence(r->d.allele[1])) continue;
            std::string raw_chr = bcf_hdr_id2name(scan.header,r->rid);
            auto chr = chromosome(raw_chr);
            if (!near_signal(windows,chr,r->pos+1,radius)) continue;
            auto k = key(chr,r->pos+1,r->d.allele[0],r->d.allele[1]);
            RefVariant v;
            v.chr = chr; v.pos = r->pos+1; v.id = r->d.id;
            v.ref = r->d.allele[0]; v.alt = r->d.allele[1]; v.key = k;
            v.row = r->rid;
            add_variant(std::move(v));
        }
    }
#if LF_WITH_PGEN
    void load_pgen_variants(PgenSource& shard, const Windows& windows, int64_t radius) {
        const auto& prefix = shard.prefix;
        std::string pvar;
        for (const auto* suffix : {".pvar", ".pvar.zst", ".pvar.gz"}) {
            if (fs::exists(prefix+suffix)) {
                require(pvar.empty(), "Ambiguous PVAR companions for " + prefix);
                pvar = prefix+suffix;
            }
        }
        require(!pvar.empty(), "Missing .pvar[.zst|.gz] for " + prefix);
        TextReader in(pvar);
        std::string line;
        bool header = false;
        while (in.getline(line)) {
            if (line.empty() || line.rfind("##", 0) == 0) continue;
            auto f = fields(line);
            if (line[0] == '#') {
                require(f.size() >= 5 && f[0] == "#CHROM" && f[1] == "POS" && f[2] == "ID" &&
                        f[3] == "REF" && f[4] == "ALT", "PVAR must begin #CHROM POS ID REF ALT: " + pvar);
                header = true; continue;
            }
            require(header && f.size() >= 5, "Malformed PVAR: " + pvar);
            require(shard.raw_count < 2147483644u, "Too many variants in PVAR");
            const uint32_t row = shard.raw_count++;
            size_t n_alleles = 2 + std::count(f[4].begin(), f[4].end(), ',');
            require(n_alleles <= 255, "PVAR exceeds pgenlib allele-count limit");
            shard.max_alleles = std::max(shard.max_alleles, static_cast<uint32_t>(n_alleles));
            if (n_alleles > 2 && shard.offsets.empty()) {
                shard.offsets.resize(static_cast<size_t>(row)+1);
                for (size_t j = 0; j < shard.offsets.size(); ++j) shard.offsets[j] = j*2;
            }
            if (!shard.offsets.empty()) shard.offsets.push_back(shard.offsets.back()+n_alleles);
            if (n_alleles != 2) { ++multiallelic_skipped; continue; }
            auto chr = chromosome(f[0]);
            size_t parsed = 0;
            int64_t pos = std::stoll(f[1], &parsed);
            require(parsed == f[1].size() && pos > 0 && pos < INT64_MAX/2, "Invalid PVAR position");
            if (!near_signal(windows, chr, pos, radius) || !sequence(f[3]) || !sequence(f[4])) continue;
            auto k = key(chr, pos, f[3], f[4]);
            RefVariant v;
            v.chr = chr; v.id = f[2]; v.ref = f[3]; v.alt = f[4]; v.key = k;
            v.pos = pos; v.row = row;
            add_variant(std::move(v));
        }
        require(header && shard.raw_count, "Empty PVAR: " + pvar);
    }
#endif
};
} // namespace lf
