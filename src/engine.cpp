#include "reference.hpp"
#include "entry.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>

namespace lf {
constexpr const char* version = "0.6.2";
constexpr size_t none = std::numeric_limits<size_t>::max();
struct Options {
    std::string input, out, regions, keep, temp_dir, missing = "skip";
    std::vector<std::string> prefixes, vcf_refs, auxiliary_inputs;
    std::string output = "-", samples = "-", expression;
    std::map<std::string,std::string> outputs;
    bool samples_file = false, exclude_filter = false, force_samples = false, record_version = true;
    std::string targets, apply_filters, index_format = "csi", command;
    bool targets_file = false;
    int regions_overlap = 1, targets_overlap = 0, compression_level = -1;
    FilterFactory filter_factory = nullptr;
    bool plugin = false, bundle_mode = false, write_index = false, regions_file = false;
    bool has_reference() const { return !prefixes.empty() || !vcf_refs.empty(); }
    double p1 = 5e-8, p2 = 0.05, lp1 = -std::log10(5e-8), indep_r2 = .6, lead_r2 = .1;
    double maf = 0, gwas_maf = 0, call_rate = 0;
    int64_t ld_bp = 1000000, merge_bp = 250000, distance_bp = 1000000;
    size_t cache_mb = 512, min_pair = 3;
    int threads = 0;
    char output_type = 'b';
    bool hardcalls = false, gwas_only = false, custom_lp = false;
};
std::string selector_path(const std::string& spec) {
    return !spec.empty() && spec.front() == '^' ? spec.substr(1) : spec;
}
void add_indexed_input(std::vector<std::string>& inputs, const std::string& path) {
    const auto delimiter = path.find(HTS_IDX_DELIM);
    if (delimiter != std::string::npos) {
        inputs.push_back(path.substr(0,delimiter));
        inputs.push_back(path.substr(delimiter+std::strlen(HTS_IDX_DELIM)));
        return;
    }
    inputs.push_back(path);
    for (const char* extension : {".csi", ".tbi"}) {
        inputs.push_back(path+extension);
        // HTSlib also looks for e.g. panel.csi beside panel.bcf.
        fs::path alternate(path);
        alternate.replace_extension(extension);
        if (fs::is_regular_file(alternate)) inputs.push_back(alternate.string());
    }
}
struct OutputFormat {
    const char* extension;
    const char* mode;
    bool indexed;
};
OutputFormat output_format(char type) {
    switch (type) {
        case 'b': return {".bcf","wb",true};
        // HTSlib maps "wbu" to level-0 BGZF; make that encoding explicit.
        case 'u': return {".bcf","wb0",false};
        case 'z': return {".vcf.gz","wz",true};
        case 'v': return {".vcf","w",false};
        default: throw std::runtime_error("output-type must be b, u, z, or v");
    }
}
std::string usage_text(bool plugin) {
    std::ostringstream out;
    out << "\nAbout: Find significant loci in GWAS-VCF summary statistics using LD clumping and FUMA-style aggregation.\n"
           "(version " << version << ")\n\n"
           "   Reads indexed GWAS-VCF/BCF summary statistics with FORMAT/LP (-log10 P),\n"
           "   treating each selected sample column as a separate phenotype.\n"
           "\n"
           "   With a reference panel, performs two-stage LD clumping to\n"
           "   identify independent significant variants and lead variants. It expands\n"
           "   LD blocks with candidate variants, joins blocks linked by LD, and merges\n"
           "   nearby blocks into genomic loci, retaining the strongest variant as each\n"
           "   locus's index. Significance, LD thresholds, and distances are configurable.\n"
           "\n"
           "   Without a reference panel, groups significant variants by physical distance.\n"
           "   Outputs selected GWAS index variants as VCF/BCF with optional locus and LD tables.\n\n";
    if (plugin)
        out << "Usage: bcftools +gwaslocusfinder [General Options] <gwas.vcf.gz> -- [Plugin Options]\n";
    else
        out << "Usage: gwaslocusfinder [options] <gwas.vcf.gz>\n"
               "       gwaslocusfinder [options] <gwas.vcf.gz> --out <prefix>\n";
    auto option = [&](const char* flag, const std::string& description) {
        out << std::left << std::setw(34) << flag << description << '\n';
    };
    out << "Options:\n";
    if (plugin) out << "   general options and the GWAS input file must precede --\n\n";
    if (plugin) {
        option("   -e, --exclude <expr>", "exclude GWAS sites for which the expression is true");
        option("   -i, --include <expr>", "select GWAS sites for which the expression is true");
    }
    option("   -f, --apply-filters <list>", "require at least one listed GWAS FILTER (e.g. \"PASS,.\")");
    option("   -r, --regions <region>", "restrict GWAS signal discovery to comma-separated regions");
    option("   -R, --regions-file <file>", "restrict GWAS signal discovery to regions listed in a file");
    option("       --regions-overlap 0|1|2", "include POS (0), record (1), or variant (2) overlaps [1]");
    option("   -t, --targets [^]<region>", "restrict GWAS signal discovery to targets; ^ excludes");
    option("   -T, --targets-file [^]<file>", "restrict GWAS signal discovery to targets in a file; ^ excludes");
    option("       --targets-overlap 0|1|2", "include POS (0), record (1), or variant (2) overlaps [0]");
    option("   -s, --samples [^]<list>", "GWAS phenotypes to include (or exclude with ^) [all]");
    option("   -S, --samples-file [^]<file>", "file of GWAS phenotypes to include (or exclude with ^)");
    option("       --force-samples", "only warn about unknown GWAS sample names");
    option("       --no-version", "do not append version and command line to the header");
    option("   -o, --output <file>", "write the set of index variants to a VCF/BCF file [standard output]");
    option("   -O, --output-type u|b|v|z[0-9]", std::string("u/b: un/compressed BCF, v/z: un/compressed VCF [") + (plugin ? "v]" : "b]"));
    out << "                                  0-9: compression level; a number alone keeps the output type\n";
    option("       --threads <int>", "use INT worker threads for HTSlib I/O [0]");
    option("   -W, --write-index[=FMT]", "index all named variant outputs: csi or tbi [off]");
    option("   -h, --help", "print this help message");
    option("   -V, --version", "print version and exit");
    if (!plugin) option("       --out <prefix>", "write the standalone TSV and variant output bundle");
    out << (plugin ? "\nPlugin options:\n" : "\nAnalysis options:\n");
    if (plugin) out << "   plugin-specific options must follow --\n\n";
#if LF_WITH_PGEN
    option("       --pfile <prefix>", "PGEN prefix or .pgen file with PVAR and PSAM; repeatable");
    option("       --pfile-list <file>", "PGEN paths, one per line, relative to the list directory");
#endif
    option("       --vcf-ref <file>", "indexed reference VCF/BCF; repeatable");
    option("       --vcf-ref-list <file>", "VCF/BCF reference paths, one per line, relative to the list");
    option("       --leads <file>", "also write second-stage LD-clumped leads; - for stdout");
    option("       --independent <file>", "also write first-stage independent variants; - for stdout");
    option("       --reports-prefix <prefix>", "also write the six TSV reports");
    out << "\nReference panel options:\n";
    option("       --keep <file>", "reference individuals: IID or FID IID, optional #header");
    option("       --maf <float>", "minimum reference MAF; monomorphic variants excluded [0]");
    option("       --min-call-rate <float>", "minimum reference called-sample fraction [0]");
    option("       --min-pair-samples <int>", "minimum complete sample pairs for LD [3]");
    option("       --hardcalls", "ignore dosages; use unphased hardcall correlation");
    option("       --missing-ref <string>", "policy for significant variants absent from panel: skip, error [skip]");
    option("       --cache-mb <int>", "decoded reference genotype cache in MiB [512]");
    option("       --temp-dir <dir>", "temporary reference storage directory [TMPDIR or .]");
    out << "\nLocus definition options:\n";
    option("       --p1 <float>", "significant when GWAS P <= threshold [5e-8]");
    option("       --log10p-threshold <float>", "significant when LP >= threshold; alternative to --p1");
    option("       --p2 <float>", "GWAS-tagged candidates require P < threshold [0.05]");
    option("       --indep-r2 <float>", "first-stage clumping and candidate r2 cutoff [0.6]");
    option("       --lead-r2 <float>", "second-stage clumping and LD locus-link cutoff [0.1]");
    option("       --ld-kb <float>", "maximum pairwise LD distance in kb [1000]");
    option("       --merge-kb <float>", "maximum gap between LD locus intervals in kb [250]");
    option("       --distance-kb <float>", "significant-variant gap without a reference panel, in kb [1000]");
    option("       --gwas-maf <float>", "per-phenotype GWAS MAF filter for signals (FORMAT/AF) [0: off]");
    option("       --gwas-only", "exclude reference-only candidates");
    out << "\nDefinitions:\n"
           "   independent: significant variants retained by the first LD-clumping stage\n"
           "                (--p1 [5e-8], --indep-r2 [0.6]).\n"
           "   leads:       a subset of independent variants retained by stricter LD clumping\n"
           "                (--lead-r2 [0.1]); several leads may belong to one locus.\n"
           "   index:       the strongest variant in each aggregated genomic locus; in LD\n"
           "                mode, index variants are a subset of leads.\n"
           "   A locus combines candidate-variant blocks linked by LD or proximity (--merge-kb).\n"
           "   Independent refers to LD pruning within --ld-kb, without conditional association testing.\n";
#if !LF_WITH_PGEN
    out << "   PGEN support is disabled in this build; use a VCF/BCF reference.\n";
#endif
    out << "\nExamples:\n";
    if (plugin) {
        out << "      bcftools +gwaslocusfinder -Ob -o index.bcf -W trait.bcf -- --vcf-ref reference.bcf\n"
               "      bcftools +gwaslocusfinder -Ob -o index.bcf -W trait.bcf -- --vcf-ref reference.bcf \\\n"
               "        --leads leads.bcf --independent independent.bcf\n"
               "      bcftools +gwaslocusfinder -Ou -o - trait.bcf -- --vcf-ref reference.bcf | bcftools view -h\n"
               "      bcftools +gwaslocusfinder -s TRAIT1,TRAIT2 -i 'FILTER=\"PASS\"' -Oz9 -o loci.vcf.gz -W=tbi \\\n"
               "        trait.bcf -- --vcf-ref reference.bcf --keep eur.keep\n"
               "      bcftools +gwaslocusfinder -Ob -o loci.bcf -W trait.bcf -- --distance-kb 1000\n";
    } else {
        out << "      gwaslocusfinder -Ob -o index.bcf -W trait.bcf --vcf-ref reference.bcf\n"
               "      gwaslocusfinder trait.bcf --out trait.loci --vcf-ref reference.bcf\n";
    }
    out << "\nNotes:\n"
           "   GWAS input and VCF/BCF references must be local, indexed BGZF VCF or BCF files.\n"
           "   All selected GWAS phenotypes are analyzed independently; --keep selects reference individuals.\n"
           "   Output is the union of selected variants, with FORMAT/LF_LOCUS per phenotype.\n"
           "   LF_START/LF_END give locus bounds; LF_INDEX_POS gives its strongest variant's position.\n"
           "   New locus annotations are written only to FORMAT, including with a single phenotype.\n"
           "   Join variant outputs and reports by phenotype and LF_LOCUS (SAMPLE and LOCUS in TSVs).\n"
           "   Index variants are always written to -o or stdout.\n"
           "   --leads/--independent add optional streams, sharing -O and -W.\n"
           "   Unselected phenotype/variant combinations have missing FORMAT values; TSVs include SAMPLE.\n"
           "   Region/target selection limits discovery; candidate expansion may extend outside it.\n"
           "   LD supports biallelic variants on chromosomes 1-22 and X, using diploid encoding.\n"
           "   Reference GT must be diploid when called; DS must use the 0-2 scale, including on X.\n"
           "   LD uses equally weighted samples, without sex-specific scaling or PAR remapping.\n";
    if (!plugin)
        out << "   Bundle mode writes TSVs plus index/lead variants; compressed outputs get CSI indexes.\n";
    out << '\n';
    return out.str();
}
double number(const std::string& s, const std::string& name) {
    size_t used = 0;
    double n = std::stod(s, &used);
    require(used == s.size() && std::isfinite(n), "Invalid number for " + name);
    return n;
}
struct HelpDone { int status = 0; };
Options parse(int argc, char** argv, bool plugin) {
    if (argc == 1) { std::cerr << usage_text(plugin); throw HelpDone{1}; }
    Options o; o.plugin = plugin;
    if (plugin) o.output_type = 'v';
    bool output_set = false, reports_set = false;
    bool p1_set = false, plugin_section = false, sample_selection_set = false;
    const std::unordered_set<std::string> general_options = {
        "-o", "--output", "-O", "--output-type", "-r", "--regions",
        "-R", "--regions-file", "--threads", "-W", "--write-index", "-V", "--version",
        "-i", "--include", "-e", "--exclude", "-s", "--samples", "-S", "--samples-file",
        "-f", "--apply-filters", "-t", "--targets", "-T", "--targets-file",
        "--regions-overlap", "--targets-overlap", "--force-samples", "--no-version"
    };
    const std::unordered_set<std::string> locus_options = {
        "--pfile", "--pfile-list", "--vcf-ref", "--vcf-ref-list", "--reports-prefix",
        "--leads", "--independent",
        "--pheno", "--p1", "--log10p-threshold", "--p2", "--indep-r2", "--lead-r2",
        "--ld-kb", "--merge-kb", "--distance-kb", "--maf", "--gwas-maf", "--min-call-rate",
        "--min-pair-samples", "--keep", "--gwas-only", "--hardcalls", "--missing-ref",
        "--cache-mb", "--temp-dir"
    };
    auto integer = [](double n, const std::string& name) -> size_t {
        require(n >= 1 && n <= 1000000000 && std::floor(n) == n, "Invalid positive integer for " + name);
        return static_cast<size_t>(n);
    };
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--") {
            require(!plugin || !plugin_section, "Only one '--' separator is allowed");
            plugin_section = true;
            continue;
        }
        if (a == "-" || a.empty() || a[0] != '-') {
            require(!plugin || !plugin_section, "GWAS INPUT must appear before '--'");
            require(o.input.empty(), "Only one input GWAS file is supported");
            o.input = a; continue;
        }
        // Normalize attached arguments before checking which side of -- owns
        // the option. The value itself (including dashes or '=') stays intact.
        std::string v;
        bool attached = false;
        const auto equal = a.find('=');
        if (a.rfind("--",0) == 0 && equal != std::string::npos) {
            v = a.substr(equal+1); a.resize(equal); attached = true;
        } else if (a.size() > 2 && a[0] == '-' && std::string("oOrRWiesSftT").find(a[1]) != std::string::npos) {
            v = a.substr(2); a.resize(2); attached = true;
            if (a == "-W" && !v.empty() && v[0] == '=') v.erase(0,1);
        }
        if (a == "--help" || a == "-h" || a == "-?") {
            require(!attached, a + " does not take a value");
            (plugin ? std::cerr : std::cout) << usage_text(plugin); throw HelpDone{};
        }
        require(a != "--emit", "--emit was removed; index output is automatic; --leads FILE and --independent FILE add optional outputs");
        require(a != "--index", "--index was removed; index output is automatic (use -o FILE or stdout)");
#if !LF_WITH_PGEN
        require(a != "--pfile" && a != "--pfile-list",
                "PGEN support is disabled; rebuild with WITH_PGEN=1 or use --vcf-ref");
#endif
        const bool general = general_options.count(a);
        require(general || locus_options.count(a) || a == "--out", "Unknown option: " + a);
        if (plugin) {
            require(a != "--out", "--out is standalone-only; use --reports-prefix after '--' for TSV reports");
            if (general) require(!plugin_section, "General option " + a + " must appear before '--'");
            else require(plugin_section, "Plugin-specific option " + a + " must appear after '--'");
        }
        if (a == "--version" || a == "-V") {
            require(!attached, a + " does not take a value");
            std::cout << version << '\n'; throw HelpDone{};
        }
        if (a == "-W" || a == "--write-index") {
            require(!attached || v == "csi" || v == "tbi", "Unsupported index format; use csi or tbi");
            o.write_index = true; o.index_format = attached ? v : "csi"; continue;
        }
        if (a == "--no-version" || a == "--force-samples") {
            require(!attached, a + " does not take a value");
            if (a == "--no-version") o.record_version = false; else o.force_samples = true;
            continue;
        }
        if (a == "--hardcalls" || a == "--gwas-only") {
            require(!attached, a + " does not take a value");
            if (a == "--hardcalls") o.hardcalls = true; else o.gwas_only = true;
            continue;
        }
        if (!attached) {
            require(i+1 < argc && std::string(argv[i+1]) != "--", "Missing value for " + a);
            v = argv[++i];
        }
        if (a == "--out") { o.out = v; o.bundle_mode = true; }
        else if (a == "-o" || a == "--output") { o.output = v; output_set = true; }
        else if (a == "--leads" || a == "--independent") {
            require(!v.empty(), a + " requires a file or '-' for stdout");
            require(o.outputs.emplace(a.substr(2),v).second, "Use " + a + " only once");
        }
        else if (a == "--reports-prefix") { o.out = v; reports_set = true; }
        else if (a == "--output-type" || a == "-O") {
            require(!v.empty(), "output-type must be u, b, v, or z, optionally followed by 0-9");
            size_t offset = 0;
            if (std::string("buzv").find(v[0]) != std::string::npos) { o.output_type = v[0]; offset = 1; }
            if (offset < v.size()) {
                require(v.size() == offset+1 && v[offset] >= '0' && v[offset] <= '9',
                        "output-type must be u, b, v, or z, optionally followed by 0-9 (or a level alone)");
                o.compression_level = v[offset]-'0';
            }
        }
        else if (a == "-i" || a == "--include" || a == "-e" || a == "--exclude") {
            require(plugin, "bcftools expressions require the plugin; pass -i/-e before '--'");
            require(o.expression.empty(), "Use only one of -i/--include or -e/--exclude");
            require(!v.empty(), "Filter expression must not be empty");
            o.expression = v; o.exclude_filter = a == "-e" || a == "--exclude";
        }
        else if (a == "-s" || a == "--samples" || a == "-S" || a == "--samples-file" || a == "--pheno") {
            require(!plugin || a != "--pheno", "Use -s NAME before '--' to select GWAS phenotypes");
            require(!sample_selection_set, "Use only one GWAS sample selection (-s or -S)");
            require(!v.empty(), "GWAS sample selection must not be empty");
            o.samples = v; o.samples_file = a == "-S" || a == "--samples-file"; sample_selection_set = true;
            if (o.samples_file) o.auxiliary_inputs.push_back(selector_path(v));
        }
        else if (a == "--regions" || a == "-r") { o.regions = v; o.regions_file = false; }
        else if (a == "--regions-file" || a == "-R") {
            o.regions = v; o.regions_file = true;
            add_indexed_input(o.auxiliary_inputs, v);
        }
        else if (a == "--targets" || a == "-t") { o.targets = v; o.targets_file = false; }
        else if (a == "--targets-file" || a == "-T") {
            o.targets = v; o.targets_file = true;
            add_indexed_input(o.auxiliary_inputs, selector_path(v));
        }
        else if (a == "--apply-filters" || a == "-f") {
            require(!v.empty(), "--apply-filters requires a nonempty FILTER list"); o.apply_filters = v;
        }
        else if (a == "--regions-overlap" || a == "--targets-overlap") {
            require(v == "0" || v == "1" || v == "2", a + " must be 0, 1, or 2");
            if (a == "--regions-overlap") o.regions_overlap = v[0]-'0'; else o.targets_overlap = v[0]-'0';
        }
        else if (a == "--keep") {
            o.keep = v;
            o.auxiliary_inputs.push_back(v);
        }
        else if (a == "--missing-ref") o.missing = v;
        else if (a == "--pfile") o.prefixes.push_back(v);
        else if (a == "--vcf-ref") o.vcf_refs.push_back(v);
        else if (a == "--pfile-list" || a == "--vcf-ref-list") {
            auto& paths = a == "--pfile-list" ? o.prefixes : o.vcf_refs;
            o.auxiliary_inputs.push_back(v);
            TextReader in(v);
            std::string line;
            size_t start = paths.size();
            while (in.getline(line)) {
                auto left = line.find_first_not_of(" \t\r");
                if (left == std::string::npos || line[left] == '#') continue;
                line = line.substr(left, line.find_last_not_of(" \t\r")-left+1);
                fs::path path(line);
                if (path.is_relative()) path = fs::absolute(v).parent_path()/path;
                paths.push_back(path.lexically_normal().string());
            }
            require(paths.size() > start, "Empty reference list: " + v);
        }
        else if (a == "--p1") { o.p1 = number(v,a); p1_set = true; }
        else if (a == "--p2") o.p2 = number(v,a);
        else if (a == "--log10p-threshold") { o.lp1 = number(v,a); o.custom_lp = true; }
        else if (a == "--indep-r2") o.indep_r2 = number(v,a);
        else if (a == "--lead-r2") o.lead_r2 = number(v,a);
        else if (a == "--maf") o.maf = number(v,a);
        else if (a == "--gwas-maf") o.gwas_maf = number(v,a);
        else if (a == "--min-call-rate") o.call_rate = number(v,a);
        else if (a == "--temp-dir") o.temp_dir = v;
        else if (a == "--cache-mb") o.cache_mb = integer(number(v,a),a);
        else if (a == "--min-pair-samples") o.min_pair = integer(number(v,a),a);
        else if (a == "--threads") {
            double n = number(v,a);
            require(n >= 0 && n <= 256 && std::floor(n) == n, "--threads must be an integer in [0,256]");
            o.threads = n;
        }
        else if (a == "--ld-kb" || a == "--merge-kb" || a == "--distance-kb") {
            double n = number(v,a)*1000;
            require(n >= 0 && n <= 1e12 && n == std::floor(n), "Distance must be a nonnegative whole number of bp");
            if (a == "--ld-kb") o.ld_bp = n;
            else if (a == "--merge-kb") o.merge_bp = n;
            else o.distance_bp = n;
        }
        else throw std::runtime_error("Unknown option: " + a);
    }
    require(!o.input.empty(), "Indexed GWAS INPUT is required (see --help)");
    require(!o.bundle_mode || (!output_set && o.outputs.empty() && !reports_set),
            "Use --out alone for a standalone bundle, or -o/explicit output destinations with --reports-prefix");
    require(!o.bundle_mode || !o.out.empty(), "--out requires a nonempty prefix");
    if (!o.bundle_mode) o.outputs.emplace("index",o.output);
    size_t stdout_count = 0, named_count = 0;
    for (const auto& output : o.outputs) {
        if (output.second == "-") ++stdout_count; else ++named_count;
        require(o.has_reference() || output.first == "index", "--" + output.first + " requires a reference panel");
    }
    require(stdout_count <= 1, "Only one output selection can use '-' for stdout");
    require(o.prefixes.empty() || o.vcf_refs.empty(), "Choose PGEN or VCF/BCF references, not both");
    require(!o.write_index || (!o.bundle_mode && named_count && output_format(o.output_type).indexed),
            "-W requires at least one named variant output with -O b or -O z");
    require(!o.write_index || o.index_format != "tbi" || o.output_type == 'z', "TBI indexing requires -Oz (BGZF VCF)");
    require(o.output_type != 'u' || o.compression_level < 0, "Compression level cannot be set on uncompressed BCF (-Ou)");
    require(o.input != "-", "Input must be a rereadable file, not stdin");
    require(!p1_set || !o.custom_lp, "Use --p1 or --log10p-threshold, not both");
    require(o.p1 > 0 && o.p1 <= 1 && o.p2 > 0 && o.p2 <= 1, "P thresholds must be in (0,1]");
    if (!o.custom_lp) o.lp1 = -std::log10(o.p1);
    require(o.lp1 >= 0 && o.lp1 >= -std::log10(o.p2), "Significance threshold must be at least as stringent as --p2");
    require(o.lead_r2 > 0 && o.lead_r2 <= o.indep_r2 && o.indep_r2 <= 1, "Require 0 < lead-r2 <= indep-r2 <= 1");
    require(o.maf >= 0 && o.maf <= .5 && o.gwas_maf >= 0 && o.gwas_maf <= .5, "MAF must be in [0,0.5]");
    require(o.call_rate >= 0 && o.call_rate <= 1, "min-call-rate must be in [0,1]");
    require(o.min_pair >= 3, "min-pair-samples must be >= 3");
    require(o.threads <= 256, "--threads must be <= 256");
    require(o.missing == "skip" || o.missing == "error", "missing-ref must be skip or error");
    require(o.has_reference() || o.keep.empty(), "--keep requires a reference panel");
    return o;
}
struct Signal {
    std::string chr, key, id;
    int64_t pos;
    int rid;
    double lp;
    Record record;
    size_t reference = none;
};
struct Exclusion { std::string key, id, reason; double lp; };
struct Association { size_t independent, reference; double r2; };
struct Independent {
    size_t signal, lead, locus = 0;
    int64_t start, end;
    bool is_lead;
};
struct Locus {
    std::string chr;
    int rid;
    int64_t start, end;
    std::vector<size_t> members;
    size_t candidates = 0;
};
struct Stats {
    size_t tagging_records = 0;
    double tagging_seconds = 0, candidate_expansion_seconds = 0;
    size_t scanned = 0, unsupported = 0, invalid_lp = 0, undefined_candidate_ld = 0, filtered = 0;
};
struct Analysis {
    std::string sample;
    int sample_index = 0;
    Stats stats;
    std::vector<Signal> signals;
    std::vector<Exclusion> excluded;
    std::vector<Independent> independent;
    std::string associations;
    std::vector<Locus> loci;
    size_t reference_phenotype = 0;
    const Signal& index_signal(const Locus& locus) const {
        // Locus members retain the significance ordering of the signals.
        return signals[independent[locus.members.front()].signal];
    }
};
std::unique_ptr<RecordFilter> make_gwas_filter(const Options& o, bcf_hdr_t* header) {
    if (o.expression.empty()) return nullptr;
    require(o.filter_factory, "bcftools filter evaluator is unavailable");
    return o.filter_factory(header,o.expression);
}
std::vector<int> gwas_filter_ids(const Options& o, bcf_hdr_t* header) {
    std::vector<int> ids;
    if (o.apply_filters.empty()) return ids;
    std::istringstream input(o.apply_filters);
    for (std::string name; std::getline(input,name,',');) {
        require(!name.empty(), "Empty FILTER name in --apply-filters");
        if (name == ".") { ids.push_back(-1); continue; }
        int id = bcf_hdr_id2int(header,BCF_DT_ID,name.c_str());
        require(id >= 0 && bcf_hdr_idinfo_exists(header,BCF_HL_FLT,id), "Unknown GWAS FILTER name: " + name);
        ids.push_back(id);
    }
    require(o.apply_filters.back() != ',', "Empty FILTER name in --apply-filters");
    return ids;
}
bool passes_filter(const Options& o, RecordFilter* filter, const std::vector<int>& ids, bcf1_t* record) {
    if (!ids.empty()) {
        bcf_unpack(record,BCF_UN_FLT);
        bool pass = !record->d.n_flt && std::find(ids.begin(),ids.end(),-1) != ids.end();
        for (int i = 0; i < record->d.n_flt; ++i)
            pass |= std::find(ids.begin(),ids.end(),record->d.flt[i]) != ids.end();
        if (!pass) return false;
    }
    return !filter || (filter->matches(record) != o.exclude_filter);
}
std::vector<int> select_gwas_samples(bcf_hdr_t* header, const Options& o) {
    std::vector<int> selected;
    if (o.samples == "-" && !o.samples_file) {
        for (int i = 0; i < bcf_hdr_nsamples(header); ++i) selected.push_back(i);
        return selected;
    }
    bool exclude = o.samples[0] == '^';
    auto spec = selector_path(o.samples);
    require(!spec.empty(), "Empty GWAS sample selection");
    int n = 0;
    char** names = hts_readlist(spec.c_str(),o.samples_file,&n);
    require(names, "Cannot read GWAS sample selection: " + spec);
    std::vector<std::string> requested;
    for (int i = 0; i < n; ++i) { requested.emplace_back(names[i]); free(names[i]); }
    free(names);
    std::set<int> seen;
    for (const auto& name : requested) {
        int i = bcf_hdr_id2int(header,BCF_DT_SAMPLE,name.c_str());
        if (i < 0 && o.force_samples) { std::cerr << "Warning: unknown GWAS sample " << name << "; skipping\n"; continue; }
        require(i >= 0, "Unknown phenotype/GWAS sample: " + name + "; use --force-samples to skip unknown names");
        require(seen.insert(i).second, "Duplicate GWAS sample selection: " + name);
        if (!exclude) selected.push_back(i);
    }
    if (exclude) for (int i = 0; i < bcf_hdr_nsamples(header); ++i)
        if (!seen.count(i)) selected.push_back(i);
    require(!selected.empty(), "No GWAS samples selected");
    return selected;
}
using Nearby = std::map<std::string, std::multimap<int64_t, size_t>>;
template<class F> void nearby(const Nearby& map, const Signal& s, int64_t radius, F f) {
    auto chr = map.find(s.chr);
    if (chr == map.end()) return;
    auto it = chr->second.lower_bound(s.pos-radius);
    for (; it != chr->second.end() && it->first <= s.pos+radius; ++it) f(it->second);
}
constexpr std::array<std::pair<const char*,const char*>,4> annotations = {{
    {"LF_LOCUS","Locus number within this GWAS phenotype; join outputs by phenotype and locus number"},
    {"LF_START","Start of the parent genomic locus on this chromosome (1-based inclusive)"},
    {"LF_END","End of the parent genomic locus on this chromosome (1-based inclusive)"},
    {"LF_INDEX_POS","Position of the strongest variant in the parent locus on this chromosome (1-based)"}
}};
void add_locus_annotations(bcf_hdr_t* header) {
    for (const auto& tag : annotations) {
        const std::string name = tag.first;
        int id = bcf_hdr_id2int(header, BCF_DT_ID, tag.first);
        if (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, id)) {
            require(bcf_hdr_id2type(header, BCF_HL_FMT, id) == BCF_HT_INT &&
                    bcf_hdr_id2length(header, BCF_HL_FMT, id) == BCF_VL_FIXED &&
                    bcf_hdr_id2number(header, BCF_HL_FMT, id) == 1,
                    "Incompatible FORMAT/" + name + ": expected Number=1,Type=Integer");
            continue;
        }
        const std::string definition = "##FORMAT=<ID=" + name + ",Number=1,Type=Integer,Description=\"" + tag.second + "\">";
        require(bcf_hdr_append(header, definition.c_str()) == 0, "Cannot add FORMAT/" + name);
    }
    require(bcf_hdr_sync(header) == 0,"Cannot sync locus annotations");
}

std::vector<Analysis> read_signals(const Options& o, Header& header, Stats& stats) {
    GwasReader in(o.input,o.regions,o.threads,o.regions_file,o.regions_overlap,o.targets,o.targets_file,o.targets_overlap);
    auto filter = make_gwas_filter(o,in.header);
    const auto filter_ids = gwas_filter_ids(o,in.header);
    const auto samples = select_gwas_samples(in.header,o);
    std::vector<Analysis> analyses(samples.size());
    std::vector<std::unordered_set<std::string>> seen(samples.size());
    for (size_t i = 0; i < samples.size(); ++i) {
        analyses[i].sample = in.header->samples[samples[i]];
        analyses[i].sample_index = samples[i];
    }
    const double threshold = static_cast<float>(o.lp1);
    std::vector<double> lp, af;
    while (auto* r = in.next()) {
        ++stats.scanned;
        // Like bcftools view, site expressions see the original input samples,
        // before selection/subsetting. No reference records enter this filter.
        if (!passes_filter(o,filter.get(),filter_ids,r)) { ++stats.filtered; continue; }
        if (r->n_allele != 2 || !sequence(r->d.allele[0]) || !sequence(r->d.allele[1])) {
            ++stats.unsupported; continue;
        }
        in.scalars(r,"LP",samples,lp);
        if (o.gwas_maf > 0) in.scalars(r,"AF",samples,af);
        auto chr = chromosome(bcf_hdr_id2name(in.header,r->rid));
        auto k = key(chr,r->pos+1,r->d.allele[0],r->d.allele[1]);
        Record copy;
        for (size_t i = 0; i < samples.size(); ++i) {
            auto& a = analyses[i];
            if (std::isnan(lp[i]) || lp[i] < 0) { ++a.stats.invalid_lp; continue; }
            if (lp[i] < threshold) continue;
            if (o.has_reference() && !autosome(chr) && chr != "X") {
                a.excluded.push_back({k,r->d.id,"unsupported_chromosome",lp[i]}); continue;
            }
            if (o.gwas_maf > 0 && (!std::isfinite(af[i]) || af[i] < 0 || af[i] > 1 || std::min(af[i],1-af[i]) < o.gwas_maf)) {
                a.excluded.push_back({k,r->d.id,"gwas_af_missing_or_below_maf",lp[i]}); continue;
            }
            require(seen[i].insert(k).second, "Duplicate significant GWAS allele pair for " + a.sample + ": " + k);
            if (!copy) copy = Record(bcf_dup(r),RecordDeleter{});
            require(copy != nullptr,"Out of memory copying GWAS record");
            a.signals.push_back({chr,k,r->d.id,r->pos+1,r->rid,lp[i],copy,none});
        }
    }
    header.reset(bcf_hdr_dup(in.header));
    require(header != nullptr,"Cannot copy GWAS header");
    add_locus_annotations(header.get());
    for (auto& a : analyses) std::sort(a.signals.begin(),a.signals.end(),[](const auto& x,const auto& y) {
        if (x.lp != y.lp) return x.lp > y.lp;
        return std::tie(x.rid,x.pos,x.key) < std::tie(y.rid,y.pos,y.key);
    });
    return analyses;
}
void tag_reference(const Options& o, const std::vector<Analysis>& analyses, Reference& ref,
                   const std::string& regions, Stats& stats) {
    // Reference keys use POS even for indels. Record overlap could visit a
    // deletion in two disjoint windows or tag a key outside both windows.
    GwasReader in(o.input,regions,o.threads,false,0);
    auto filter = make_gwas_filter(o,in.header);
    const auto filter_ids = gwas_filter_ids(o,in.header);
    std::vector<int> samples;
    for (auto& a : analyses) samples.push_back(a.sample_index);
    std::vector<double> lp;
    while (auto* r = in.next()) {
        ++stats.tagging_records;
        const bool eligible = r->n_allele == 2 && passes_filter(o,filter.get(),filter_ids,r);
        bool decoded = false;
        // Filtered GWAS variants remain tagged as present, with missing LP,
        // so they cannot re-enter as reference-only candidates.
        for (int alt = 1; alt < r->n_allele; ++alt) {
            if (!sequence(r->d.allele[0]) || !sequence(r->d.allele[alt])) continue;
            auto chr = chromosome(bcf_hdr_id2name(in.header,r->rid));
            auto k = key(chr,r->pos+1,r->d.allele[0],r->d.allele[alt]);
            auto found = ref.store.find(chr,r->pos+1,k);
            if (!found) continue;
            ref.store.tag(*found,k);
            if (eligible) {
                if (!decoded) { in.scalars(r,"LP",samples,lp); decoded = true; }
                ref.store.set_lp(*found,lp);
            }
        }
    }
}
DisjointSet find_ld(const Options& o, Analysis& analysis, Reference& ref) {
    auto& signals = analysis.signals;
    auto& independent = analysis.independent;
    std::unique_ptr<BinaryFile> associations;
    if (!o.out.empty()) {
        analysis.associations = ref.scratch.file();
        associations = std::make_unique<BinaryFile>(analysis.associations,"wb");
    }
    const float candidate_threshold = static_cast<float>(-std::log10(o.p2));
    require(ref.subset.size() >= o.min_pair, "Reference subset smaller than --min-pair-samples");
    // First-stage clumping retains independent significant variants.
    Nearby chosen;
    for (size_t i = 0; i < signals.size(); ++i) {
        auto& s = signals[i];
        auto hit = ref.store.find(s.chr,s.pos,s.key);
        if (!hit) {
            require(o.missing != "error", "Significant variant missing from reference: " + s.key);
            analysis.excluded.push_back({s.key,s.id,"missing_reference",s.lp}); continue;
        }
        s.reference = *hit;
        if (!ref.passes(s.reference)) {
            analysis.excluded.push_back({s.key,s.id,"reference_maf_callrate_or_variance",s.lp}); continue;
        }
        auto g = ref.genotype(s.reference);
        bool clumped = false;
        nearby(chosen, s, o.ld_bp, [&](size_t j) {
            if (clumped) return;
            auto h = ref.genotype(signals[independent[j].signal].reference);
            double r2 = r_squared(*g,*h,o.min_pair);
            require(std::isfinite(r2), "Undefined LD between significant variants (missingness/zero pairwise variance): " +
                    s.key + " and " + signals[independent[j].signal].key);
            if (r2 >= o.indep_r2) clumped = true;
        });
        if (!clumped) {
            size_t j = independent.size();
            independent.push_back({i,j,0,s.pos,s.pos,true});
            chosen[s.chr].emplace(s.pos,j);
        }
    }
    // Second-stage clumping selects leads and links their candidate blocks.
    DisjointSet components(independent.size());
    Nearby previous;
    for (size_t i = 0; i < independent.size(); ++i) {
        auto& a = independent[i];
        auto& s = signals[a.signal];
        auto g = ref.genotype(s.reference);
        size_t owner = none;
        nearby(previous, s, o.ld_bp, [&](size_t j) {
            auto h = ref.genotype(signals[independent[j].signal].reference);
            double r2 = r_squared(*g,*h,o.min_pair);
            require(std::isfinite(r2), "Undefined LD between independent variants: " + s.key);
            if (r2 >= o.lead_r2) {
                components.join(i,j);
                if (independent[j].is_lead) owner = std::min(owner,j);
            }
        });
        if (owner != none) { a.is_lead = false; a.lead = owner; }
        previous[s.chr].emplace(s.pos,i);
        const auto candidate_start = std::chrono::steady_clock::now();
        ref.store.neighborhood(s.chr, s.pos, o.ld_bp, [&](const RefVariant& v) {
            size_t j = v.index;
            if (j != s.reference) {
                const bool tagged = ref.store.tagged(j);
                if (tagged && !(ref.store.lp(j,analysis.reference_phenotype) > candidate_threshold)) return;
                if (!tagged && o.gwas_only) return;
                if (!ref.passes(j)) return;
            }
            double r2 = j == s.reference ? 1 : r_squared(*g,*ref.genotype(j),o.min_pair);
            if (!std::isfinite(r2)) { ++analysis.stats.undefined_candidate_ld; return; }
            if (r2 >= o.indep_r2) {
                if (associations) associations->write(Association{i,j,r2});
                a.start = std::min(a.start,v.pos); a.end = std::max(a.end,v.pos);
            }
        });
        analysis.stats.candidate_expansion_seconds +=
            std::chrono::duration<double>(std::chrono::steady_clock::now()-candidate_start).count();
    }
    // Union all low-threshold LD edges first, then merge component envelopes.
    // This preserves chains A--B--C even when A and C are both lead variants.
    std::map<size_t,Interval> spans;
    for (size_t i = 0; i < independent.size(); ++i) {
        auto& a = independent[i];
        size_t root = components.find(i);
        auto found = spans.find(root);
        if (found == spans.end()) spans.emplace(root,Interval{signals[a.signal].chr,a.start,a.end,i});
        else { found->second.start = std::min(found->second.start,a.start); found->second.end = std::max(found->second.end,a.end); }
    }
    std::vector<Interval> intervals;
    for (const auto& item : spans) intervals.push_back(item.second);
    merge_intervals(std::move(intervals),o.merge_bp,components);
    if (associations) associations->close();
    return components;
}
DisjointSet find_distance(const Options& o, Analysis& analysis) {
    DisjointSet components(analysis.signals.size());
    std::vector<Interval> intervals;
    for (size_t i = 0; i < analysis.signals.size(); ++i) {
        const auto& signal = analysis.signals[i];
        analysis.independent.push_back({i,i,0,signal.pos,signal.pos,true});
        intervals.push_back({signal.chr,signal.pos,signal.pos,i});
    }
    merge_intervals(std::move(intervals),o.distance_bp,components);
    return components;
}
std::vector<Locus> collect_loci(const std::vector<Signal>& signals, std::vector<Independent>& independent,
                                DisjointSet& dsu) {
    std::map<size_t,Locus> grouped;
    for (size_t i = 0; i < independent.size(); ++i) {
        const auto& a = independent[i];
        const auto& s = signals[a.signal];
        auto root = dsu.find(i);
        auto it = grouped.find(root);
        if (it == grouped.end()) it = grouped.emplace(root,Locus{s.chr,s.rid,a.start,a.end,{},{}}).first;
        auto& loc = it->second;
        loc.start = std::min(loc.start,a.start); loc.end = std::max(loc.end,a.end);
        loc.members.push_back(i);
    }
    std::vector<Locus> out;
    for (auto& kv : grouped) out.push_back(std::move(kv.second));
    std::sort(out.begin(),out.end(),[](const auto& a,const auto& b) {
        return std::tie(a.rid,a.start,a.end) < std::tie(b.rid,b.start,b.end);
    });
    for (size_t i = 0; i < out.size(); ++i)
        for (size_t j : out[i].members) independent[j].locus = i+1;
    return out;
}
// Candidate edges are only retained when reports were requested. Sorting
// (final locus, reference ID) counts unique variants without large locus sets.
struct LocusCandidate {
    size_t locus = 0, reference = 0;
    bool operator<(const LocusCandidate& b) const { return std::tie(locus,reference) < std::tie(b.locus,b.reference); }
    bool read(BinaryFile& file) { return file.read(*this); }
    void write(BinaryFile& file) const { file.write(*this); }
    size_t memory_bytes() const { return sizeof(*this); }
};
void count_candidates(Analysis& analysis, ScratchDirectory& scratch) {
    if (analysis.associations.empty()) return;
    ExternalSort<LocusCandidate> sorter(scratch);
    BinaryFile associations(analysis.associations,"rb");
    Association edge{};
    while (associations.read(edge))
        sorter.add({analysis.independent[edge.independent].locus,edge.reference});
    std::optional<LocusCandidate> previous;
    sorter.finish([&](const LocusCandidate& value) {
        if (!previous || *previous < value) ++analysis.loci[value.locus-1].candidates;
        previous = value;
    });
}
// Output files and both possible index sidecars must have distinct paths.
// Resolve aliases so one output cannot replace another output or an input.
bool same_path(const std::string& a, const std::string& b) {
    if (fs::weakly_canonical(a) == fs::weakly_canonical(b)) return true;
    return fs::exists(a) && fs::exists(b) && fs::equivalent(a,b);
}
struct OutputPaths {
    std::vector<std::string> inputs, outputs;
    explicit OutputPaths(const Options& o) : inputs(o.auxiliary_inputs) {
        add_indexed_input(inputs, o.input);
        for (const auto& path : o.vcf_refs) {
            add_indexed_input(inputs, path);
            // Reference opens the canonical filename, including its sidecars.
            if (fs::is_symlink(path)) add_indexed_input(inputs, fs::canonical(path).string());
        }
        for (auto path : o.prefixes) {
            if (ends_with(path,".pgen")) path.resize(path.size()-5);
            for (const auto* suffix : {".pgen",".pgen.pgi",".psam",".pvar",".pvar.zst",".pvar.gz"}) inputs.push_back(path+suffix);
        }
    }
    void add(const std::string& path, bool allow_stream = false) {
        require(!path.empty(), "Output path must not be empty");
        for (const auto& input : inputs)
            require(!same_path(path,input), "Output/index path collides with input: " + path);
        for (const auto& output : outputs)
            require(!same_path(path,output), "Output/index path collides with another output or report: " + path);
        require(!fs::exists(path) || fs::is_regular_file(path)
                || (allow_stream && (fs::is_fifo(path) || fs::is_character_file(path))),
                "Output is not a regular file or supported stream: " + path);
        require(fs::is_directory(fs::absolute(path).parent_path()), "Output directory does not exist: " + path);
        outputs.push_back(path);
    }
    void variant(const std::string& path, bool indexed = false) {
        if (path == "-") return;
        // HTSlib strips ##idx## before opening the data file, including on writes.
        require(path.find(HTS_IDX_DELIM) == std::string::npos,
                "Explicit ##idx## filenames are unsupported for variant outputs; use -W to write an index");
        add(path,!indexed); add(path+".csi"); add(path+".tbi");
    }
};
struct ReportOutputs {
    std::string prefix;
    ReportOutputs(const std::string& p, char output_type, bool reports_only, OutputPaths& paths) : prefix(p) {
        for (const auto* suffix : {".loci.tsv",".leads.tsv",".independent.tsv",".candidates.tsv",".excluded.tsv",".run.tsv"})
            paths.add(prefix+suffix);
        if (!reports_only) for (const auto* selection : {".index",".leads"})
            paths.variant(prefix+selection+output_format(output_type).extension,output_format(output_type).indexed);
    }
    std::string path(const std::string& suffix) const { return prefix+suffix; }
    std::ofstream stream(const std::string& suffix) const {
        std::ofstream out(path(suffix));
        out.exceptions(std::ios::badbit|std::ios::failbit);
        out << std::setprecision(10);
        return out;
    }
    // Invalidate the previous completion marker before any new output is written.
    void begin_write() const { fs::remove(path(".run.tsv")); }
};
struct VariantOutput {
    std::string path, emit;
    bool indexed;
    VariantOutput(const std::string& destination, bool index, const std::string& selection, OutputPaths& paths)
        : path(destination), emit(selection), indexed(index && destination != "-") {
        paths.variant(destination,indexed);
    }
};
// Mask every FORMAT field for unselected sample/variant cells. Missing
// vectors are encoded as one missing value followed by vector-end sentinels.
void mask_unselected(bcf_hdr_t* header, bcf1_t* record, const std::vector<int32_t>& loci) {
    if (std::find(loci.begin(),loci.end(),bcf_int32_missing) == loci.end()) return;
    bcf_unpack(record,BCF_UN_FMT);
    std::vector<std::pair<std::string,int>> tags;
    for (int i = 0; i < record->n_fmt; ++i) {
        int id = record->d.fmt[i].id;
        tags.emplace_back(bcf_hdr_int2id(header,BCF_DT_ID,id),bcf_hdr_id2type(header,BCF_HL_FMT,id));
    }
    struct Buffers {
        int32_t* integer = nullptr; float* real = nullptr; char* text = nullptr;
        int ni = 0, nf = 0, nc = 0;
        ~Buffers() { free(integer); free(real); free(text); }
    } buf;
    for (const auto& tag : tags) {
        int n = 0, status = -1;
        if (tag.first == "GT" || tag.second == BCF_HT_INT) {
            n = tag.first == "GT" ? bcf_get_genotypes(header,record,&buf.integer,&buf.ni)
                : bcf_get_format_int32(header,record,tag.first.c_str(),&buf.integer,&buf.ni);
            require(n > 0 && n % loci.size() == 0,"Invalid integer FORMAT field: " + tag.first);
            size_t stride = n/loci.size();
            for (size_t i = 0; i < loci.size(); ++i) if (loci[i] == bcf_int32_missing) {
                buf.integer[i*stride] = tag.first == "GT" ? bcf_gt_missing : bcf_int32_missing;
                std::fill(buf.integer+i*stride+1,buf.integer+(i+1)*stride,bcf_int32_vector_end);
            }
            status = tag.first == "GT" ? bcf_update_genotypes(header,record,buf.integer,n)
                : bcf_update_format_int32(header,record,tag.first.c_str(),buf.integer,n);
        } else if (tag.second == BCF_HT_REAL) {
            n = bcf_get_format_float(header,record,tag.first.c_str(),&buf.real,&buf.nf);
            require(n > 0 && n % loci.size() == 0,"Invalid float FORMAT field: " + tag.first);
            size_t stride = n/loci.size();
            for (size_t i = 0; i < loci.size(); ++i) if (loci[i] == bcf_int32_missing) {
                bcf_float_set_missing(buf.real[i*stride]);
                for (size_t j = 1; j < stride; ++j) bcf_float_set_vector_end(buf.real[i*stride+j]);
            }
            status = bcf_update_format_float(header,record,tag.first.c_str(),buf.real,n);
        } else if (tag.second == BCF_HT_STR) {
            n = bcf_get_format_char(header,record,tag.first.c_str(),&buf.text,&buf.nc);
            require(n > 0 && n % loci.size() == 0,"Invalid string FORMAT field: " + tag.first);
            size_t stride = n/loci.size();
            for (size_t i = 0; i < loci.size(); ++i) if (loci[i] == bcf_int32_missing) {
                std::fill(buf.text+i*stride,buf.text+(i+1)*stride,'\0');
                buf.text[i*stride] = '.';
            }
            status = bcf_update_format_char(header,record,tag.first.c_str(),buf.text,n);
        }
        require(status == 0,"Cannot mask FORMAT/" + tag.first);
    }
}
void write_selected_variants(const std::string& path, bcf_hdr_t* header,
                             const std::vector<Analysis>& analyses, const std::string& emit,
                             const Options& o, bool index_output) {
    std::vector<char*> samples;
    for (const auto& a : analyses) samples.push_back(const_cast<char*>(a.sample.c_str()));
    std::vector<int> imap(samples.size());
    Header h(bcf_hdr_subset(header,samples.size(),samples.data(),imap.data()));
    require(h && bcf_hdr_nsamples(h.get()) == static_cast<int>(samples.size()),"Cannot subset output GWAS header");
    if (o.record_version) {
        std::string meta = "##gwaslocusfinderVersion=" + std::string(version);
        require(bcf_hdr_append(h.get(),meta.c_str()) == 0,"Cannot annotate output header");
        meta = "##gwaslocusfinderCommand=" + o.command;
        require(bcf_hdr_append(h.get(),meta.c_str()) == 0,"Cannot annotate output command");
    }
    require(bcf_hdr_sync(h.get()) == 0,"Cannot sync output header");
    struct Selected { const Signal* signal; std::vector<int32_t> loci; };
    std::vector<Selected> selection;
    std::unordered_map<bcf1_t*,size_t> seen;
    auto add = [&](size_t sample, const Signal& signal, size_t locus) {
        require(locus <= INT32_MAX,"Too many loci for FORMAT/LF_LOCUS");
        auto inserted = seen.emplace(signal.record.get(),selection.size());
        if (inserted.second) selection.push_back({&signal,std::vector<int32_t>(analyses.size(),bcf_int32_missing)});
        selection[inserted.first->second].loci[sample] = locus;
    };
    for (size_t i = 0; i < analyses.size(); ++i) {
        const auto& a = analyses[i];
        if (emit == "index") {
            for (size_t j = 0; j < a.loci.size(); ++j) add(i,a.index_signal(a.loci[j]),j+1);
        } else if (o.has_reference()) for (const auto& b : a.independent)
            if (emit == "independent" || b.is_lead) add(i,a.signals[b.signal],b.locus);
    }
    std::stable_sort(selection.begin(),selection.end(),[](const auto& a,const auto& b) {
        return std::tie(a.signal->rid,a.signal->pos,a.signal->key) < std::tie(b.signal->rid,b.signal->pos,b.signal->key);
    });
    const auto format = output_format(o.output_type);
    std::string mode = format.mode;
    if (o.compression_level >= 0 && (o.output_type == 'b' || o.output_type == 'z'))
        mode += static_cast<char>('0'+o.compression_level);
    htsFile* out = bcf_open(path.c_str(),mode.c_str());
    require(out,"Cannot create output VCF/BCF: " + path);
    try {
        // The destination is open for direct writing. Its previous indexes no
        // longer describe it, even if a later write or indexing step fails.
        if (path != "-") for (const auto* suffix : {".csi",".tbi"}) fs::remove(path+suffix);
        if (o.threads > 0) require(hts_set_threads(out,o.threads) == 0,"Cannot set output threads");
        require(bcf_hdr_write(out,h.get()) == 0,"Cannot write output header");
        for (const auto& item : selection) {
            Record r(bcf_dup(item.signal->record.get()),RecordDeleter{});
            require(r != nullptr,"Cannot copy selected GWAS record");
            require(bcf_subset(header,r.get(),samples.size(),imap.data()) == 0,"Cannot subset GWAS record");
            mask_unselected(h.get(),r.get(),item.loci);
            std::array<std::vector<int32_t>,4> values;
            values[0] = item.loci;
            for (size_t j = 1; j < values.size(); ++j) values[j].assign(analyses.size(),bcf_int32_missing);
            for (size_t i = 0; i < analyses.size(); ++i) if (item.loci[i] != bcf_int32_missing) {
                const auto& analysis = analyses[i];
                const auto& locus = analysis.loci[item.loci[i]-1];
                const auto& index = analysis.index_signal(locus);
                require(locus.start <= INT32_MAX && locus.end <= INT32_MAX && index.pos <= INT32_MAX,
                        "Locus coordinates exceed the VCF FORMAT Integer range");
                values[1][i] = locus.start; values[2][i] = locus.end; values[3][i] = index.pos;
            }
            for (size_t j = 0; j < annotations.size(); ++j) {
                const auto* tag = annotations[j].first;
                require(bcf_update_format_int32(h.get(),r.get(),tag,values[j].data(),values[j].size()) == 0,
                        "Cannot annotate FORMAT/" + std::string(tag));
            }
            require(bcf_write(out,h.get(),r.get()) == 0,"Cannot write selected record");
        }
    } catch (...) { bcf_close(out); throw; }
    require(bcf_close(out) == 0,"Cannot close output VCF/BCF: " + path);
    if (index_output && format.indexed)
        require(bcf_index_build3(path.c_str(),nullptr,o.index_format == "tbi" ? 0 : 14,o.threads) == 0,"Cannot index output VCF/BCF: " + path);
}
void write_outputs(ReportOutputs& reports, const Options& o, bcf_hdr_t* header,
                   const std::vector<Analysis>& analyses, const Stats& stats,
                   Reference* ref, double elapsed, const std::string& command) {
    const std::string method = o.has_reference() ? "ld" : "distance";
    size_t n_signals = 0, n_excluded = 0, n_loci = 0, n_independent = 0;
    std::string phenotypes;
    for (const auto& a : analyses) {
        n_signals += a.signals.size(); n_excluded += a.excluded.size();
        n_loci += a.loci.size(); n_independent += a.independent.size();
        if (!phenotypes.empty()) phenotypes += ',';
        phenotypes += a.sample;
    }
    {
        auto out = reports.stream(".loci.tsv");
        out << "SAMPLE\tLOCUS\tCHROM\tSTART\tEND\tINDEX_KEY\tINDEX_ID\tINDEX_POS\tINDEX_LP\tN_INDEPENDENT\tN_LEAD\tN_CANDIDATE\tMETHOD\n";
        for (const auto& a : analyses) for (size_t i = 0; i < a.loci.size(); ++i) {
            const auto& loc = a.loci[i];
            const auto& s = a.index_signal(loc);
            size_t leads = 0;
            for (auto j : loc.members) leads += a.independent[j].is_lead;
            out << a.sample << '\t' << i+1 << '\t' << loc.chr << '\t' << loc.start << '\t' << loc.end << '\t'
                << s.key << '\t' << s.id << '\t' << s.pos << '\t' << s.lp << '\t';
            if (ref) out << loc.members.size() << '\t' << leads << '\t' << loc.candidates;
            else out << ".\t.\t.";
            out << '\t' << method << '\n';
        }
        out.close();
    }
    for (const auto& suffix : {".independent.tsv",".leads.tsv"}) {
        auto out = reports.stream(suffix);
        out << "SAMPLE\tLOCUS\tKEY\tID\tCHROM\tPOS\tLP\tBLOCK_START\tBLOCK_END\tLEAD_KEY\tLOCUS_START\tLOCUS_END\tINDEX_KEY\tINDEX_ID\tINDEX_POS\tMETHOD\n";
        if (ref) for (const auto& a : analyses) for (const auto& b : a.independent) {
            if (std::string(suffix) == ".leads.tsv" && !b.is_lead) continue;
            const auto& s = a.signals[b.signal];
            out << a.sample << '\t' << b.locus << '\t' << s.key << '\t' << s.id << '\t' << s.chr << '\t' << s.pos << '\t'
                << s.lp << '\t' << b.start << '\t' << b.end << '\t'
                << a.signals[a.independent[b.lead].signal].key << '\t';
            const auto& locus = a.loci[b.locus-1];
            const auto& index = a.index_signal(locus);
            out << locus.start << '\t' << locus.end << '\t' << index.key << '\t' << index.id << '\t' << index.pos << '\t' << method << '\n';
        }
        out.close();
    }
    {
        auto out = reports.stream(".candidates.tsv");
        out << "SAMPLE\tLOCUS\tINDEPENDENT_KEY\tKEY\tID\tCHROM\tPOS\tREF\tALT\tR2\tMAF\tCALL_RATE\tGWAS_LP\tSOURCE\n";
        if (ref) for (const auto& a : analyses) {
            if (a.associations.empty()) continue;
            BinaryFile associations(a.associations,"rb");
            Association b{};
            while (associations.read(b)) {
                const auto v = ref->store.get(b.reference);
                auto qc = ref->store.qc(b.reference);
                require(qc.has_value(),"Missing candidate QC");
                out << a.sample << '\t' << a.independent[b.independent].locus << '\t' << a.signals[a.independent[b.independent].signal].key << '\t'
                    << v.key << '\t' << v.id << '\t' << v.chr << '\t' << v.pos << '\t' << v.ref << '\t' << v.alt << '\t'
                    << b.r2 << '\t' << qc->maf << '\t' << qc->call_rate << '\t';
                float lp = ref->store.lp(b.reference,a.reference_phenotype);
                if (std::isnan(lp)) out << '.'; else out << lp;
                out << '\t' << (ref->store.tagged(b.reference) ? "gwas" : "reference_only") << '\n';
            }
        }
        out.close();
    }
    {
        auto out = reports.stream(".excluded.tsv");
        out << "SAMPLE\tKEY\tID\tLP\tREASON\n";
        for (const auto& a : analyses) for (const auto& v : a.excluded)
            out << a.sample << '\t' << v.key << '\t' << v.id << '\t' << v.lp << '\t' << v.reason << '\n';
        out.close();
    }
    if (o.bundle_mode) {
        const auto format = output_format(o.output_type);
        for (const char* emit : {"index","leads"})
            write_selected_variants(reports.path(std::string(".")+emit+format.extension),header,analyses,emit,o,true);
    }
    {
        auto out = reports.stream(".run.tsv");
        out << "PARAMETER\tVALUE\nversion\t" << version << "\nmethod\t" << method
            << "\noutput_type\t" << o.output_type << "\ncompression_level\t" << o.compression_level
            << "\nindex_format\t" << o.index_format << "\nrecord_version\t" << o.record_version
            << "\noutput_mode\t" << (o.bundle_mode ? "bundle" : "streams")
            << "\ninput\t" << o.input << "\nphenotype\t" << phenotypes << "\ngwas_samples\t" << analyses.size()
            << "\nsamples_selection\t" << o.samples << "\nsamples_file\t" << o.samples_file
            << "\ngwas_filter_mode\t" << (o.expression.empty() ? "none" : o.exclude_filter ? "exclude" : "include")
            << "\ngwas_filter_expression\t" << o.expression << "\ngwas_apply_filters\t" << o.apply_filters
            << "\nforce_samples\t" << o.force_samples << "\nrecords_filtered\t" << stats.filtered
            << "\nregions\t" << o.regions << "\nregions_file\t" << o.regions_file
            << "\nregions_overlap\t" << o.regions_overlap << "\ntargets\t" << o.targets
            << "\ntargets_file\t" << o.targets_file << "\ntargets_overlap\t" << o.targets_overlap
            << "\nlog10p_threshold\t" << o.lp1 << "\np2\t" << o.p2
            << "\nindep_r2\t" << o.indep_r2 << "\nlead_r2\t" << o.lead_r2 << "\nld_bp\t" << o.ld_bp
            << "\nmerge_bp\t" << o.merge_bp << "\ndistance_bp\t" << o.distance_bp
            << "\nreference_maf\t" << o.maf << "\ngwas_maf\t" << o.gwas_maf << "\nmin_call_rate\t" << o.call_rate
            << "\nmin_pair_samples\t" << o.min_pair << "\nkeep\t" << o.keep
            << "\nld_measure\t" << (o.hardcalls ? "unphased_hardcall_r2" : "unphased_dosage_r2")
            << "\ninclude_reference_only\t" << !o.gwas_only << "\nmissing_reference_policy\t" << o.missing
            << "\nreference_storage\t" << (ref ? "temporary_binary" : ".") << "\ntemp_dir\t" << o.temp_dir
            << "\ncache_mb\t" << o.cache_mb << "\nhts_threads\t" << o.threads << "\nhtslib_version\t" << hts_version()
            << "\npgen_support\t" << LF_WITH_PGEN
            << "\npgenlib_revision\t" << (LF_WITH_PGEN ? "4544f92e374e784a6ffd97fa88239962b81f2b62" : ".")
            << "\nrecords_scanned\t" << stats.scanned << "\nunsupported_gwas_records\t" << stats.unsupported
            << "\ninvalid_or_missing_lp\t" << stats.invalid_lp << "\nsignificant_records\t" << n_signals
            << "\nexcluded_signals\t" << n_excluded << "\nloci\t" << n_loci
            << "\nundefined_candidate_ld\t" << stats.undefined_candidate_ld << "\nelapsed_seconds\t" << elapsed << '\n';
        if (ref) {
            out << "reference_backend\t" << (o.vcf_refs.empty() ? "pgen" : "vcf") << "\nreference_samples\t" << ref->subset.size()
                << "\ntagging_records\t" << stats.tagging_records << "\ntagging_seconds\t" << stats.tagging_seconds
                << "\ncandidate_expansion_seconds\t" << stats.candidate_expansion_seconds
                << "\nreference_variants_in_windows\t" << ref->store.size() << "\nreference_duplicates_skipped\t" << ref->duplicates_skipped
                << "\nreference_variants_qc_evaluated\t" << ref->qc_counts.evaluated
                << "\nreference_variants_qc_passed\t" << ref->qc_counts.passed
                << "\nreference_variants_qc_filtered\t" << ref->qc_counts.filtered()
                << "\nreference_variants_qc_not_evaluated\t" << ref->qc_not_evaluated()
                << "\nreference_multiallelic_skipped\t" << ref->multiallelic_skipped << "\ngenotype_reads\t" << ref->genotype_reads
                << "\nindependent_variants\t" << n_independent << '\n';
            for (const auto& shard : ref->shards) out << "reference_prefix\t" << shard->prefix << '\n';
        }
        if (o.bundle_mode) {
            for (const auto* selection : {"index","leads"})
                out << selection << "_output\t" << reports.path(std::string(".")+selection+output_format(o.output_type).extension) << '\n';
        } else for (const auto& output : o.outputs) out << output.first << "_output\t" << output.second << '\n';
        out << "command\t" << command << "\nstatus\tcomplete\n";
        out.close();
    }
}
int run(int argc, char** argv, bool plugin, FilterApi filters) {
    const auto start = std::chrono::steady_clock::now();
    Options o = parse(argc,argv,plugin);
    o.filter_factory = filters.compile;
    if (!o.expression.empty()) {
        auto paths = filters.input_files(o.expression);
        o.auxiliary_inputs.insert(o.auxiliary_inputs.end(),paths.begin(),paths.end());
    }
    std::ostringstream command;
    for (int i = 0; i < argc; ++i) { if (i) command << ' '; command << std::quoted(argv[i]); }
    for (char c : command.str()) {
        if (c == '\n') o.command += "\\n";
        else if (c == '\r') o.command += "\\r";
        else if (c == '\t') o.command += "\\t";
        else o.command += c;
    }
    OutputPaths paths(o);
    std::unique_ptr<ReportOutputs> reports;
    if (!o.out.empty()) reports = std::make_unique<ReportOutputs>(o.out,o.output_type,!o.bundle_mode,paths);
    std::vector<VariantOutput> outputs;
    for (const auto& output : o.outputs) outputs.emplace_back(output.second,o.write_index,output.first,paths);
    Header header;
    Stats stats;
    std::unique_ptr<Reference> ref;
    auto analyses = read_signals(o,header,stats);
    std::cerr << "Read " << stats.scanned << " records; analyzing " << analyses.size() << " GWAS phenotypes\n";
    Windows windows;
    for (const auto& a : analyses) for (const auto& s : a.signals) windows[s.chr].push_back(s.pos);
    for (auto& item : windows) {
        std::sort(item.second.begin(),item.second.end());
        item.second.erase(std::unique(item.second.begin(),item.second.end()),item.second.end());
    }
    if (o.has_reference() && (!windows.empty() || !o.vcf_refs.empty())) {
        ref = std::make_unique<Reference>(o.vcf_refs.empty() ? o.prefixes : o.vcf_refs,windows,o.ld_bp,o.keep,o.cache_mb,o.hardcalls,!o.vcf_refs.empty(),o.threads,o.temp_dir,analyses.size(),ReferenceFilters{o.maf,o.call_rate});
        if (ref->duplicates_skipped)
            std::cerr << "Warning: Skipped " << ref->duplicates_skipped
                      << " duplicate reference rows in queried windows; kept the first occurrence of each allele-pair key"
                      << " (showing at most 10 duplicate details).\n";
        std::cerr << "Reference: " << ref->store.size() << " nearby variants before reference QC, " << ref->subset.size() << " selected samples\n";
        const auto tagging_start = std::chrono::steady_clock::now();
        const auto regions = signal_regions(header.get(), windows, o.ld_bp, &ref->store);
        if (!regions.empty()) tag_reference(o,analyses,*ref,regions,stats);
        ref->store.finish_tagging();
        stats.tagging_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-tagging_start).count();
    }
    for (size_t phenotype = 0; phenotype < analyses.size(); ++phenotype) {
        auto& a = analyses[phenotype];
        a.reference_phenotype = phenotype;
        auto components = ref ? find_ld(o,a,*ref) : find_distance(o,a);
        a.loci = collect_loci(a.signals,a.independent,components);
        if (ref) count_candidates(a,ref->scratch);
        stats.candidate_expansion_seconds += a.stats.candidate_expansion_seconds;
        stats.invalid_lp += a.stats.invalid_lp;
        stats.undefined_candidate_ld += a.stats.undefined_candidate_ld;
        std::cerr << a.sample << ": " << a.loci.size() << " loci; " << a.excluded.size() << " excluded significant variants\n";
    }
    if (ref) {
        const auto& qc = ref->qc_counts;
        std::cerr << "Reference QC: " << qc.evaluated << " variants evaluated, " << qc.passed
                  << " passed, " << qc.filtered() << " filtered out; " << ref->qc_not_evaluated() << " not evaluated (genotypes not needed)\n";
    }
    double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if (reports) reports->begin_write();
    for (const auto& output : outputs)
        write_selected_variants(output.path,header.get(),analyses,output.emit,o,output.indexed);
    if (reports) write_outputs(*reports,o,header.get(),analyses,stats,ref.get(),elapsed,o.command);
    if (stats.unsupported) std::cerr << "Skipped " << stats.unsupported << " multiallelic/symbolic GWAS records; normalize/split before analysis\n";
    return 0;
}
int entry(int argc, char** argv, bool plugin, FilterApi filters) {
    try { return run(argc,argv,plugin,filters); }
    catch (const HelpDone& e) { return e.status; }
    catch (const std::exception& e) { std::cerr << "gwaslocusfinder: " << e.what() << '\n'; return 1; }
}
} // namespace lf
