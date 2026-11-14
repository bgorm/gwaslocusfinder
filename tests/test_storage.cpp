#include "storage.hpp"
#include <iostream>
#include <map>
#include <random>

using namespace lf;
static RefVariant variant(size_t i) {
    RefVariant v;
    v.chr = i < 6000 ? "1" : "2";
    v.pos = i < 1200 ? 100 : static_cast<int64_t>(i);
    v.ref = "A"; v.alt = "A";
    for (size_t n = i, j = 0; j < 8; ++j, n >>= 2) v.alt += "ACGT"[n & 3];
    v.id = "v"+std::to_string(i); v.key = key(v.chr,v.pos,v.ref,v.alt);
    v.shard = i%3; v.row = i+123; v.index = i;
    return v;
}
static void check(bool value) { require(value,"Storage regression check failed"); }
int main() {
    try {
        ScratchDirectory parent("");
        const std::string directory = parent.file();
        fs::create_directory(directory);
        {
            ScratchDirectory scratch(directory);
            // Hundreds of tiny runs force multiple merge passes and duplicate
            // ties. Compare the complete result against an in-memory oracle.
            ExternalSort<RefVariant,RefOrder> sorter(scratch,2048);
            std::vector<RefVariant> expected;
            for (size_t i = 0; i < 2400; ++i) {
                auto v = variant((i*37)%1200); v.index = i;
                expected.push_back(v); sorter.add(v);
            }
            std::sort(expected.begin(),expected.end(),RefOrder{});
            size_t seen = 0;
            sorter.finish([&](const RefVariant& v) {
                const auto& e = expected.at(seen++);
                check(v.key == e.key && v.index == e.index && v.row == e.row && v.shard == e.shard);
            });
            check(seen == expected.size());

            ReferenceStore store(scratch,3);
            std::vector<size_t> order(10000);
            std::iota(order.begin(),order.end(),0);
            std::mt19937 rng(42); std::shuffle(order.begin(),order.end(),rng);
            for (auto i : order) store.add(variant(i));
            for (size_t i = 0; i < 100; ++i) { auto v = variant(i); v.row = 999999; store.add(v); }
            size_t duplicates = 0;
            store.build([&](const RefVariant& first,const RefVariant& later) {
                check(first.row != 999999 && later.row == 999999); ++duplicates;
            });
            check(duplicates == 100 && store.size() == 10000);
            // High blocks first: intervening annotation blocks are file holes.
            for (auto it = order.rbegin(); it != order.rend(); ++it) {
                auto v = variant(*it); auto id = store.find(v.chr,v.pos,v.key);
                check(id.has_value());
                const auto actual = store.get(*id);
                check(actual.row == v.row && actual.shard == v.shard && actual.id == v.id);
                if (*it%3 == 0) {
                    store.tag(*id,v.key);
                    if (*it%2 == 0) store.set_lp(*id,{0,std::numeric_limits<double>::quiet_NaN(),static_cast<double>(*it)});
                }
                if (*it%5 == 0) store.set_qc(*id,{.25,.75,static_cast<double>(*it)});
            }
            store.finish_tagging();
            // Random revisits exceed the eight-block cache many times. Missing
            // LP must remain distinct from zero, including filtered GWAS tags.
            for (auto i : order) {
                auto v = variant(i); size_t id = *store.find(v.chr,v.pos,v.key);
                check(store.tagged(id) == (i%3 == 0));
                if (i%6 == 0) { check(store.lp(id,0) == 0 && store.lp(id,2) == i); }
                else check(std::isnan(store.lp(id,0)) && std::isnan(store.lp(id,2)));
                check(std::isnan(store.lp(id,1)));
                auto qc = store.qc(id);
                check(qc.has_value() == (i%5 == 0));
                if (qc) check(qc->maf == .25 && qc->call_rate == .75 && qc->variance == i);
            }
            for (auto [chr,pos,radius] : std::vector<std::tuple<std::string,int64_t,int64_t>>{
                     {"1",100,0},{"1",5900,200},{"2",6000,1},{"2",9999,1},{"3",100,100}}) {
                std::vector<std::string> actual, wanted;
                store.neighborhood(chr,pos,radius,[&](const RefVariant& v) {
                    const auto key_before = v.key;
                    for (size_t i = 0; i < store.size(); i += 512) store.get(i);
                    check(v.key == key_before); // active window block remains valid
                    actual.push_back(v.key);
                });
                for (size_t i = 0; i < 10000; ++i) {
                    auto v = variant(i);
                    if (v.chr == chr && v.pos >= pos-radius && v.pos <= pos+radius) wanted.push_back(v.key);
                }
                std::sort(actual.begin(),actual.end()); std::sort(wanted.begin(),wanted.end());
                check(actual == wanted);
            }
            check(!store.find("1",1,"1:1:A:G") && !store.find("3",100,"3:100:A:G"));
            auto v = variant(0); auto id = *store.find(v.chr,v.pos,v.key);
            bool duplicate_failed = false;
            try { store.tag(id,v.key); } catch (const std::runtime_error&) { duplicate_failed = true; }
            check(duplicate_failed);
        }
        check(fs::is_empty(directory));
        // Cleanup also runs on exceptions, without touching pre-existing files.
        auto sentinel = fs::path(directory)/"keep.txt";
        { std::ofstream out(sentinel); out << "keep"; }
        try {
            ScratchDirectory scratch(directory);
            BinaryFile out(scratch.file(),"wb"); out.write(uint64_t{7});
            throw std::runtime_error("expected");
        } catch (const std::runtime_error&) {}
        check(std::distance(fs::directory_iterator(directory),fs::directory_iterator{}) == 1 && fs::exists(sentinel));
        if (fs::exists("/dev/full")) {
            bool failed = false;
            try { BinaryFile out("/dev/full","wb"); out.write(uint64_t{7}); out.flush(); }
            catch (const std::runtime_error&) { failed = true; }
            check(failed);
        }
        std::cout << "Temporary storage checks passed\n";
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
