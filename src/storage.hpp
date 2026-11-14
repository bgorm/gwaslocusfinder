#pragma once
#include "io.hpp"
#include <cerrno>
#include <cstring>
#include <cstdlib>
#include <list>
#include <optional>
#include <queue>
#include <type_traits>
#include <unistd.h>

namespace lf {
// All scratch files are private to this run. These are native binary files,
// not a persistent interchange format; no database or extra library is needed.
class ScratchDirectory {
    fs::path path_;
    size_t serial_ = 0;
public:
    explicit ScratchDirectory(std::string directory) {
        if (directory.empty()) {
            const char* env = std::getenv("TMPDIR");
            directory = env && *env ? env : ".";
        }
        std::string pattern = (fs::absolute(directory)/"gwaslocusfinder.XXXXXX").string();
        auto* created = mkdtemp(pattern.data());
        require(created, "Cannot create temporary directory in " + directory + ": " + std::strerror(errno));
        path_ = pattern;
    }
    ScratchDirectory(const ScratchDirectory&) = delete;
    ~ScratchDirectory() { std::error_code ec; fs::remove_all(path_,ec); }
    std::string file() { return (path_/(std::to_string(serial_++)+".bin")).string(); }
};
class BinaryFile {
    FILE* file_ = nullptr;
    std::string path_;
public:
    BinaryFile(const std::string& path, const char* mode) : path_(path) {
        file_ = fopen(path.c_str(),mode);
        require(file_, "Cannot open temporary file " + path + ": " + std::strerror(errno));
    }
    BinaryFile(const BinaryFile&) = delete;
    ~BinaryFile() { if (file_) fclose(file_); }
    void flush() { require(fflush(file_) == 0,"Cannot flush temporary file " + path_ + " (check free disk space)"); }
    void close() {
        FILE* f = file_; file_ = nullptr;
        require(fclose(f) == 0,"Cannot close temporary file " + path_);
    }
    void seek(uint64_t offset) {
        require(offset <= INT64_MAX && fseeko(file_,static_cast<off_t>(offset),SEEK_SET) == 0,
                "Cannot seek temporary file " + path_);
    }
    void resize(uint64_t bytes) {
        flush();
        require(bytes <= INT64_MAX && ftruncate(fileno(file_),static_cast<off_t>(bytes)) == 0,
                "Cannot size temporary file " + path_);
    }
    uint64_t tell() {
        auto offset = ftello(file_);
        require(offset >= 0,"Cannot locate position in temporary file " + path_);
        return offset;
    }
    bool read_bytes(void* data, size_t bytes, bool allow_eof = false) {
        if (!bytes) return true;
        size_t n = fread(data,1,bytes,file_);
        require(!ferror(file_) && (n == bytes || (!n && allow_eof && feof(file_))),
                "Cannot read temporary file " + path_);
        return n != 0;
    }
    void write_bytes(const void* data, size_t bytes) {
        require(fwrite(data,1,bytes,file_) == bytes,"Cannot write temporary file " + path_ + " (check free disk space)");
    }
    template<class T> bool read(T& value) {
        static_assert(std::is_trivially_copyable_v<T>);
        return read_bytes(&value,sizeof(value),true);
    }
    template<class T> void write(const T& value) {
        static_assert(std::is_trivially_copyable_v<T>);
        write_bytes(&value,sizeof(value));
    }
    void read_string(std::string& value) {
        uint64_t n = 0;
        require(read(n) && n <= value.max_size(),"Invalid temporary string in " + path_);
        value.resize(n);
        read_bytes(value.data(),value.size());
    }
    void write_string(const std::string& value) {
        write(static_cast<uint64_t>(value.size())); write_bytes(value.data(),value.size());
    }
};
struct RefVariant {
    std::string chr, id, ref, alt, key;
    int64_t pos = 0;
    size_t index = 0, shard = 0;
    uint32_t row = 0;
    bool read(BinaryFile& file) {
        if (!file.read(index)) return false;
        require(file.read(pos) && file.read(shard) && file.read(row),"Truncated temporary reference record");
        file.read_string(chr); file.read_string(id); file.read_string(ref); file.read_string(alt);
        key = lf::key(chr,pos,ref,alt);
        return true;
    }
    void write(BinaryFile& file) const {
        file.write(index); file.write(pos); file.write(shard); file.write(row);
        file.write_string(chr); file.write_string(id); file.write_string(ref); file.write_string(alt);
    }
    size_t memory_bytes() const { return sizeof(*this)+chr.capacity()+id.capacity()+ref.capacity()+alt.capacity()+key.capacity(); }
};
inline auto coordinate(const RefVariant& v) { return std::tie(v.chr,v.pos,v.key); }
struct RefOrder {
    bool operator()(const RefVariant& a, const RefVariant& b) const {
        if (coordinate(a) != coordinate(b)) return coordinate(a) < coordinate(b);
        return a.index < b.index; // original input order breaks duplicate ties
    }
};
// Bounded-memory external merge sort. A bounded fan-in also caps open files.
// The callback consumes the final merge directly, without another full copy.
template<class T, class Less = std::less<T>> class ExternalSort {
    ScratchDirectory& scratch_;
    std::vector<T> buffer_;
    std::vector<std::string> runs_;
    size_t bytes_ = 0, limit_;
    Less less_;
    void spill() {
        if (buffer_.empty()) return;
        std::sort(buffer_.begin(),buffer_.end(),less_);
        auto path = scratch_.file();
        BinaryFile out(path,"wb");
        for (const auto& value : buffer_) value.write(out);
        out.close(); runs_.push_back(path);
        std::vector<T>().swap(buffer_); bytes_ = 0;
    }
    template<class F> void merge(const std::vector<std::string>& paths, F visit) {
        std::vector<std::unique_ptr<BinaryFile>> files;
        std::vector<T> heads(paths.size());
        auto later = [&](size_t a, size_t b) { return less_(heads[b],heads[a]); };
        std::priority_queue<size_t,std::vector<size_t>,decltype(later)> pending(later);
        for (size_t i = 0; i < paths.size(); ++i) {
            files.push_back(std::make_unique<BinaryFile>(paths[i],"rb"));
            if (heads[i].read(*files[i])) pending.push(i);
        }
        while (!pending.empty()) {
            size_t i = pending.top(); pending.pop();
            visit(heads[i]);
            if (heads[i].read(*files[i])) pending.push(i);
        }
    }
public:
    explicit ExternalSort(ScratchDirectory& scratch, size_t bytes = 32*1024*1024)
        : scratch_(scratch), limit_(bytes) {}
    void add(T value) {
        bytes_ += value.memory_bytes(); buffer_.push_back(std::move(value));
        if (bytes_ >= limit_) spill();
    }
    template<class F> void finish(F visit) {
        if (runs_.empty()) {
            std::sort(buffer_.begin(),buffer_.end(),less_);
            for (const auto& v : buffer_) visit(v);
            std::vector<T>().swap(buffer_); bytes_ = 0;
            return;
        }
        spill();
        constexpr size_t fan_in = 32;
        while (runs_.size() > fan_in) {
            std::vector<std::string> next;
            for (size_t i = 0; i < runs_.size(); i += fan_in) {
                std::vector<std::string> group(runs_.begin()+i,runs_.begin()+std::min(runs_.size(),i+fan_in));
                auto path = scratch_.file();
                BinaryFile out(path,"wb");
                merge(group,[&](const T& value) { value.write(out); });
                out.close(); next.push_back(path);
                for (const auto& input : group) fs::remove(input);
            }
            runs_ = std::move(next);
        }
        merge(runs_,visit);
        for (const auto& input : runs_) fs::remove(input);
        runs_.clear();
    }
};
struct RefQC { double maf = 0, call_rate = 0, variance = 0; };
// Sparse coordinate index + a small LRU of metadata/annotation blocks. Variant
// IDs are sequential ordinals in the sorted file; PGEN rows stay separate.
class ReferenceStore {
    static constexpr size_t block_size = 512, cache_blocks = 8;
    static constexpr uint64_t tagged_flag = 1, lp_flag = 2, qc_flag = 4;
    struct State { uint64_t flags = 0; RefQC qc; };
    struct Index { std::string chr, key; int64_t pos; uint64_t offset; };
    struct Block {
        std::vector<RefVariant> variants;
        std::array<State,block_size> states{};
        std::vector<float> lp;
        bool state_dirty = false, lp_dirty = false;
        size_t phenotype = SIZE_MAX, memory = 0;
        std::list<size_t>::iterator position;
    };
    ExternalSort<RefVariant,RefOrder> sorter_;
    BinaryFile metadata_, states_, lp_;
    std::vector<Index> index_;
    std::list<size_t> lru_;
    std::unordered_map<size_t,std::shared_ptr<Block>> cache_;
    size_t count_ = 0, input_count_ = 0, phenotypes_, cache_bytes_ = 0;
    uint64_t lp_block_bytes_;
    bool tagging_ = true;
    static uint64_t offset(size_t block, uint64_t bytes) {
        require(bytes && block <= INT64_MAX/bytes,"Temporary reference annotation file is too large");
        return block*bytes;
    }
    void save(size_t id, Block& block) {
        if (block.state_dirty) {
            states_.seek(offset(id,sizeof(block.states)));
            states_.write_bytes(block.states.data(),sizeof(block.states));
            block.state_dirty = false;
        }
        if (block.lp_dirty) {
            lp_.seek(offset(id,lp_block_bytes_));
            lp_.write_bytes(block.lp.data(),block.lp.size()*sizeof(float));
            block.lp_dirty = false;
        }
    }
    void trim() {
        constexpr size_t byte_limit = 32*1024*1024;
        // A single oversized block (or a block pinned by a callback) may exceed
        // the cache target. All other cached blocks are evicted first.
        while (cache_.size() > 1 && (cache_.size() > cache_blocks || cache_bytes_ > byte_limit)) {
            size_t old = lru_.back(); auto b = cache_.at(old); save(old,*b);
            cache_bytes_ -= b->memory; cache_.erase(old); lru_.pop_back();
        }
    }
    void load_lp(size_t id, Block& block, size_t phenotype) {
        const size_t selected = tagging_ ? SIZE_MAX : phenotype;
        if (!block.lp.empty() && block.phenotype == selected) return;
        size_t n = tagging_ ? block_size*phenotypes_ : block_size;
        block.lp.resize(n);
        lp_.seek(offset(id,lp_block_bytes_)+(tagging_ ? 0 : phenotype*block_size*sizeof(float)));
        lp_.read_bytes(block.lp.data(),n*sizeof(float));
        block.phenotype = selected;
        // Capacity only grows within a phase; finish_tagging drops the cache.
        size_t memory = sizeof(Block)+block.lp.capacity()*sizeof(float)+block.variants.capacity()*sizeof(RefVariant);
        for (const auto& v : block.variants) memory += v.memory_bytes()-sizeof(RefVariant);
        cache_bytes_ += memory-block.memory; block.memory = memory;
        trim();
    }
    std::shared_ptr<Block> block(size_t id) {
        auto found = cache_.find(id);
        if (found != cache_.end()) {
            lru_.splice(lru_.begin(),lru_,found->second->position);
            return found->second;
        }
        auto value = std::make_shared<Block>();
        metadata_.seek(index_.at(id).offset);
        size_t n = std::min(block_size,count_-id*block_size);
        value->variants.resize(n);
        for (auto& v : value->variants) require(v.read(metadata_),"Truncated temporary reference metadata");
        states_.seek(offset(id,sizeof(value->states)));
        states_.read_bytes(value->states.data(),sizeof(value->states));
        value->memory = sizeof(Block)+value->variants.capacity()*sizeof(RefVariant);
        for (const auto& v : value->variants) value->memory += v.memory_bytes()-sizeof(RefVariant);
        lru_.push_front(id); value->position = lru_.begin(); cache_[id] = value;
        cache_bytes_ += value->memory; trim();
        return value;
    }
    size_t locate(const std::string& chr, int64_t pos, const std::string& key) const {
        auto it = std::upper_bound(index_.begin(),index_.end(),std::tie(chr,pos,key),
                                  [](const auto& value,const Index& i) { return value < std::tie(i.chr,i.pos,i.key); });
        return it == index_.begin() ? 0 : static_cast<size_t>(it-index_.begin()-1);
    }
public:
    // Scratch ownership is supplied by Reference and outlives this store.
    ReferenceStore(ScratchDirectory& scratch, size_t phenotypes)
        : sorter_(scratch), metadata_(scratch.file(),"w+b"), states_(scratch.file(),"w+b"), lp_(scratch.file(),"w+b"),
          phenotypes_(phenotypes), lp_block_bytes_(0) {
        static_assert(sizeof(State) == 32);
        require(phenotypes > 0 && phenotypes <= INT64_MAX/(block_size*sizeof(float)),"Too many phenotypes for reference annotations");
        lp_block_bytes_ = block_size*sizeof(float)*phenotypes;
    }
    void add(RefVariant value) { value.index = input_count_++; sorter_.add(std::move(value)); }
    template<class F> void build(F duplicate) {
        std::optional<RefVariant> previous;
        sorter_.finish([&](const RefVariant& input) {
            if (previous && previous->key == input.key) { duplicate(*previous,input); return; }
            previous = input;
            auto value = input; value.index = count_;
            if (count_%block_size == 0) index_.push_back({value.chr,value.key,value.pos,metadata_.tell()});
            value.write(metadata_); ++count_;
        });
        metadata_.flush();
        // Logical zero-filled files make untouched blocks cheap while ensuring
        // unexpected EOF is an error, never mistaken for missing annotations.
        states_.resize(offset(index_.size(),block_size*sizeof(State)));
        lp_.resize(offset(index_.size(),lp_block_bytes_));
    }
    void finish_tagging() {
        for (auto& entry : cache_) save(entry.first,*entry.second);
        states_.flush(); lp_.flush();
        cache_.clear(); lru_.clear(); cache_bytes_ = 0; tagging_ = false;
    }
    size_t size() const { return count_; }
    RefVariant get(size_t id) { require(id < count_,"Invalid reference ID"); return block(id/block_size)->variants[id%block_size]; }
    std::optional<size_t> find(const std::string& chr, int64_t pos, const std::string& key) {
        if (index_.empty()) return {};
        auto values = block(locate(chr,pos,key));
        auto it = std::lower_bound(values->variants.begin(),values->variants.end(),std::tie(chr,pos,key),
                                   [](const RefVariant& v,const auto& wanted) { return coordinate(v) < wanted; });
        if (it == values->variants.end() || it->key != key) return {};
        return it->index;
    }
    void tag(size_t id, const std::string& key) {
        auto b = block(id/block_size); auto& state = b->states[id%block_size];
        require(!(state.flags & tagged_flag),"Duplicate GWAS allele pair in reference window: " + key);
        state.flags |= tagged_flag; b->state_dirty = true;
    }
    void set_lp(size_t id, const std::vector<double>& lp) {
        require(tagging_ && lp.size() == phenotypes_,"Inconsistent reference phenotype tagging");
        auto b = block(id/block_size); load_lp(id/block_size,*b,0);
        for (size_t i = 0; i < phenotypes_; ++i) b->lp[i*block_size+id%block_size] = lp[i];
        b->states[id%block_size].flags |= lp_flag; b->state_dirty = true; b->lp_dirty = true;
    }
    float lp(size_t id, size_t phenotype) {
        require(phenotype < phenotypes_,"Invalid reference phenotype index");
        auto b = block(id/block_size);
        if (!(b->states[id%block_size].flags & lp_flag)) return std::numeric_limits<float>::quiet_NaN();
        load_lp(id/block_size,*b,phenotype);
        return b->lp[(tagging_ ? phenotype*block_size : 0)+id%block_size];
    }
    bool tagged(size_t id) { return block(id/block_size)->states[id%block_size].flags & tagged_flag; }
    std::optional<RefQC> qc(size_t id) {
        auto b = block(id/block_size); const auto& s = b->states[id%block_size];
        return s.flags & qc_flag ? std::optional<RefQC>(s.qc) : std::nullopt;
    }
    void set_qc(size_t id, RefQC value) {
        auto b = block(id/block_size); auto& s = b->states[id%block_size];
        s.qc = value; s.flags |= qc_flag; b->state_dirty = true;
    }
    template<class F> void neighborhood(const std::string& chr, int64_t pos, int64_t radius, F visit) {
        if (index_.empty()) return;
        for (size_t i = locate(chr,pos-radius,""); i < index_.size(); ++i) {
            auto b = block(i); // pinned while callbacks may access other blocks
            for (const auto& v : b->variants) {
                if (v.chr < chr || (v.chr == chr && v.pos < pos-radius)) continue;
                if (v.chr != chr || v.pos > pos+radius) return;
                visit(v);
            }
        }
    }
};
} // namespace lf
