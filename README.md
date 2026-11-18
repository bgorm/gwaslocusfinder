# gwaslocusfinder

Find significant loci in [GWAS-VCF](https://github.com/MRCIEU/gwas-vcf-specification)
summary statistics using LD clumping and FUMA-style aggregation.
Gwaslocusfinder identifies independent significant variants, LD-clumped lead
variants, and one index variant per aggregated genomic locus. Each GWAS sample
column is analyzed as a separate phenotype.

**gwaslocusfinder is a bcftools plugin** designed to work with Giulio Genovese's
[bcftools score suite](https://github.com/freeseek/score). It reads LD reference
panels directly from PLINK PGEN or indexed VCF/BCF,
including panels split across multiple files. Without a reference panel, it
finds loci by physical distance. Index variants are written as VCF/BCF;
lead variants, independent variants, and tabular reports are optional outputs.

- [Usage](#usage)
- [Installation](#installation)
- [Find significant loci](#find-significant-loci)
- [Worked example using Alzheimer’s disease summary statistics](#alzheimers-disease-on-chromosome-19)
- [Reference panels](#reference-panels)
- [Locus definition](#locus-definition)
- [Outputs](#outputs)
- [Clumping and thresholding for polygenic scores](#clumping-and-thresholding-for-polygenic-scores)
- [Performance and limitations](#performance-and-limitations)
- [License](#license)

## Usage

```sh
bcftools +gwaslocusfinder --help
```

<!-- BEGIN GWASLOCUSFINDER HELP -->
```text
About: Find significant loci in GWAS-VCF summary statistics using LD clumping and FUMA-style aggregation.
(version 0.6.1)

   Reads indexed GWAS-VCF/BCF summary statistics with FORMAT/LP (-log10 P),
   treating each selected sample column as a separate phenotype.

   With a reference panel, performs two-stage LD clumping to
   identify independent significant variants and lead variants. It expands
   LD blocks with candidate variants, joins blocks linked by LD, and merges
   nearby blocks into genomic loci, retaining the strongest variant as each
   locus's index. Significance, LD thresholds, and distances are configurable.

   Without a reference panel, groups significant variants by physical distance.
   Outputs selected GWAS index variants as VCF/BCF with optional locus and LD tables.

Usage: bcftools +gwaslocusfinder [General Options] <gwas.vcf.gz> -- [Plugin Options]
Options:
   general options and the GWAS input file must precede --

   -e, --exclude <expr>           exclude GWAS sites for which the expression is true
   -i, --include <expr>           select GWAS sites for which the expression is true
   -f, --apply-filters <list>     require at least one listed GWAS FILTER (e.g. "PASS,.")
   -r, --regions <region>         restrict GWAS signal discovery to comma-separated regions
   -R, --regions-file <file>      restrict GWAS signal discovery to regions listed in a file
       --regions-overlap 0|1|2    include POS (0), record (1), or variant (2) overlaps [1]
   -t, --targets [^]<region>      restrict GWAS signal discovery to targets; ^ excludes
   -T, --targets-file [^]<file>   restrict GWAS signal discovery to targets in a file; ^ excludes
       --targets-overlap 0|1|2    include POS (0), record (1), or variant (2) overlaps [0]
   -s, --samples [^]<list>        GWAS phenotypes to include (or exclude with ^) [all]
   -S, --samples-file [^]<file>   file of GWAS phenotypes to include (or exclude with ^)
       --force-samples            only warn about unknown GWAS sample names
       --no-version               do not append version and command line to the header
   -o, --output <file>            write the set of index variants to a VCF/BCF file [standard output]
   -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF [v]
                                  0-9: compression level; a number alone keeps the output type
       --threads <int>            use INT worker threads for HTSlib I/O [0]
   -W, --write-index[=FMT]        index all named variant outputs: csi or tbi [off]
   -h, --help                     print this help message
   -V, --version                  print version and exit

Plugin options:
   plugin-specific options must follow --

       --pfile <prefix>           PGEN prefix or .pgen file with PVAR and PSAM; repeatable
       --pfile-list <file>        PGEN paths, one per line, relative to the list directory
       --vcf-ref <file>           indexed reference VCF/BCF; repeatable
       --vcf-ref-list <file>      VCF/BCF reference paths, one per line, relative to the list
       --leads <file>             also write second-stage LD-clumped leads; - for stdout
       --independent <file>       also write first-stage independent variants; - for stdout
       --reports-prefix <prefix>  also write the six TSV reports

Reference panel options:
       --keep <file>              reference individuals: IID or FID IID, optional #header
       --maf <float>              minimum reference MAF; monomorphic variants excluded [0]
       --min-call-rate <float>    minimum reference called-sample fraction [0]
       --min-pair-samples <int>   minimum complete sample pairs for LD [3]
       --hardcalls                ignore dosages; use unphased hardcall correlation
       --missing-ref <string>     policy for significant variants absent from panel: skip, error [skip]
       --cache-mb <int>           decoded reference genotype cache in MiB [512]
       --temp-dir <dir>           temporary reference storage directory [TMPDIR or .]

Locus definition options:
       --p1 <float>               significant when GWAS P <= threshold [5e-8]
       --log10p-threshold <float> significant when LP >= threshold; alternative to --p1
       --p2 <float>               GWAS-tagged candidates require P < threshold [0.05]
       --indep-r2 <float>         first-stage clumping and candidate r2 cutoff [0.6]
       --lead-r2 <float>          second-stage clumping and LD locus-link cutoff [0.1]
       --ld-kb <float>            maximum pairwise LD distance in kb [1000]
       --merge-kb <float>         maximum gap between LD locus intervals in kb [250]
       --distance-kb <float>      significant-variant gap without a reference panel, in kb [1000]
       --gwas-maf <float>         per-phenotype GWAS MAF filter for signals (FORMAT/AF) [0: off]
       --gwas-only                exclude reference-only candidates

Definitions:
   independent: significant variants retained by the first LD-clumping stage
                (--p1 [5e-8], --indep-r2 [0.6]).
   leads:       a subset of independent variants retained by stricter LD clumping
                (--lead-r2 [0.1]); several leads may belong to one locus.
   index:       the strongest variant in each aggregated genomic locus; in LD
                mode, index variants are a subset of leads.
   A locus combines candidate-variant blocks linked by LD or proximity (--merge-kb).
   Independent refers to LD pruning within --ld-kb, without conditional association testing.

Examples:
      bcftools +gwaslocusfinder -Ob -o index.bcf -W trait.bcf -- --vcf-ref reference.bcf
      bcftools +gwaslocusfinder -Ob -o index.bcf -W trait.bcf -- --vcf-ref reference.bcf \
        --leads leads.bcf --independent independent.bcf
      bcftools +gwaslocusfinder -Ou -o - trait.bcf -- --vcf-ref reference.bcf | bcftools view -h
      bcftools +gwaslocusfinder -s TRAIT1,TRAIT2 -i 'FILTER="PASS"' -Oz9 -o loci.vcf.gz -W=tbi \
        trait.bcf -- --vcf-ref reference.bcf --keep eur.keep
      bcftools +gwaslocusfinder -Ob -o loci.bcf -W trait.bcf -- --distance-kb 1000

Notes:
   GWAS input and VCF/BCF references must be local, indexed BGZF VCF or BCF files.
   All selected GWAS phenotypes are analyzed independently; --keep selects reference individuals.
   Output is the union of selected variants, with FORMAT/LF_LOCUS per phenotype.
   LF_START/LF_END give locus bounds; LF_INDEX_POS gives its strongest variant's position.
   New locus annotations are written only to FORMAT, including with a single phenotype.
   Join variant outputs and reports by phenotype and LF_LOCUS (SAMPLE and LOCUS in TSVs).
   Index variants are always written to -o or stdout.
   --leads/--independent add optional streams, sharing -O and -W.
   Unselected phenotype/variant combinations have missing FORMAT values; TSVs include SAMPLE.
   Region/target selection limits discovery; candidate expansion may extend outside it.
   LD supports biallelic variants on chromosomes 1-22 and X, using diploid encoding.
   Reference GT must be diploid when called; DS must use the 0-2 scale, including on X.
   LD uses equally weighted samples, without sex-specific scaling or PAR remapping.
```
<!-- END GWASLOCUSFINDER HELP -->

## Installation

The plugin requires a C++17 compiler, GNU make, bcftools, matching HTSlib
headers (1.25 or newer), and the zlib and zstd development libraries.
GCC and Clang are supported compiler choices. On macOS, use the Xcode Command
Line Tools or another C++17 toolchain; dependencies can come from a package
manager or source builds. bcftools and the plugin must use the same architecture
(for example, both arm64 on Apple Silicon).

The bootstrap also uses Bash, curl, tar, and either `sha256sum` or `shasum`.
It downloads checksum-verified sources and builds only `gwaslocusfinder.so`.
It does not install anything globally or compile a separate copy of HTSlib.
The plugin resolves HTSlib symbols from its bcftools host.

### Default build with PGEN support

Clone the repository and build the plugin:

```sh
git clone https://github.com/bgorm/gwaslocusfinder.git
cd gwaslocusfinder
bash scripts/bootstrap.sh
export BCFTOOLS_PLUGINS="$PWD/build:${BCFTOOLS_PLUGINS:-}"
bcftools +gwaslocusfinder --help
```

This fetches HTSlib 1.25 headers and a pinned pgenlib revision, then builds
`build/gwaslocusfinder.so` with PGEN and VCF/BCF reference support. The PLINK
executable and Python are not runtime dependencies. The plugin search path
retains existing entries and bcftools' default directories when previously
unset. A full path also works: `bcftools +/path/to/gwaslocusfinder.so ...`.

The default headers match bcftools/HTSlib 1.25. For another installation,
point to matching HTSlib sources. Both version and architecture must match;
the plugin rejects a host with different bcftools or HTSlib versions.

```sh
PLUGIN_HTSLIB=/path/to/matching/htslib-source \
  bash scripts/bootstrap.sh
```

### Minimal build without pgenlib

For VCF/BCF references and distance mode only:

```sh
bash scripts/bootstrap.sh --without-pgen
```

This does **not download, compile, or link pgenlib**. It keeps the same LD
clumping, sample selection, filters, output annotations, and locus aggregation
for VCF/BCF references. `--pfile` and `--pfile-list` are unavailable and produce
a build-capability error if requested. Use `--vcf-ref` or `--vcf-ref-list`.
The help output above shows the default PGEN-enabled build.

With dependencies already available, the equivalent make command is:

```sh
make -j4 WITH_PGEN=0
```

Keep `WITH_PGEN=0` on subsequent make commands such as `make test` and
`make install`. The default is `WITH_PGEN=1`. Backend configurations use
separate object directories, so changing this flag does not reuse incompatible
objects. zlib, zstd, and matching HTSlib headers remain required in both builds.

### Compiler and dependency paths

`make` builds the plugin by default. Standard `CPPFLAGS`, `CXXFLAGS`, `LDFLAGS`,
and `LDLIBS` are honored. `pkg-config`, when available, supplies zlib/zstd flags;
otherwise the compiler's normal include and library search paths are used.
`PGENLIB` can point to a PLINK source tree's `2.0/include` directory, including
the adjacent `2.0/simde` headers needed on ARM.

For dependencies installed outside the compiler's search paths, supply their
include and library directories, for example:

```sh
make -j4 CXX=clang++ \
  CPPFLAGS="-I/path/to/dependencies/include" \
  LDFLAGS="-L/path/to/dependencies/lib" \
  PLUGIN_HTSLIB=/path/to/matching/htslib-source
```

Installed HTSlib headers can also be used. Set `PLUGIN_HTSLIB` to the directory
containing `htslib/hts.h`, and provide its version explicitly because installed
headers do not include HTSlib's `version.sh`. With an appropriate pkg-config
search path:

```sh
make -j4 WITH_PGEN=0 \
  PLUGIN_HTSLIB="$(pkg-config --variable=includedir htslib)" \
  PLUGIN_HTSLIB_VERSION="$(pkg-config --modversion htslib)"
```

Use `BCFTOOLS=/path/to/bcftools` to select a particular host. Run `make clean`
before changing compiler flags, dependencies, or host versions. On macOS the
Makefile uses a native bundle with host-symbol lookup and keeps the `.so`
extension expected by bcftools. CI builds and tests the plugin on Linux and
macOS, with and without PGEN support.

### Install and test

To install into a chosen bcftools plugin directory:

```sh
make install plugindir=/path/to/bcftools/plugins
```

`PREFIX` (default `/usr/local`) and `DESTDIR` are supported. Alternatively,
load the plugin directly from `build/` with `BCFTOOLS_PLUGINS` as above.
For a minimal installation, add `WITH_PGEN=0` to the install command.

The integration tests require Python 3 and bcftools. PGEN-enabled builds also
read a small synthetic PGEN and compare its results with the matching VCF.
The fixture is included in the repository.

```sh
make test
# For a build without pgenlib:
make test WITH_PGEN=0
```

## Find significant loci

### Prepare the GWAS input

To convert tabular GWAS summary statistics (such as in GWAS-SSF format) to
GWAS-VCF, use
[`bcftools +munge`](https://github.com/freeseek/score#convert-summary-statistics):

```sh
set -euo pipefail
bcftools +munge -c SSF \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -s TRAIT -Ou sumstats.tsv.gz | \
  bcftools view -f PASS,. -m2 -M2 -Ou | \
  bcftools sort -Ob -o trait.bcf -W
```

This example assumes GRCh38 coordinates and uses the
[GRCh38 no-alt reference used by freeseek/score](https://github.com/freeseek/score#installation),
with its `.fai` index. `-c SSF` selects GWAS-SSF column mappings; for other
tables use another `-c` preset or `-C colheaders.tsv`. `TRAIT` names the
phenotype and `trait.bcf` is the sorted, indexed output.

GWAS input must be a **local, coordinate-sorted, indexed BCF or BGZF-compressed
VCF**. A matching CSI or TBI index is required even when no regions are selected.
Plain VCF, ordinary gzip, stdin, and remote URLs are not supported as input.
For an existing GWAS-VCF, prepare an indexed BCF with:

```sh
bcftools sort -Ob -o trait.bcf trait.vcf.gz
bcftools index trait.bcf
```

Each phenotype needs `FORMAT/LP`, the negative base-10 logarithm of the P-value,
with `Type=Float` and `Number=A` or `Number=1`. Missing or invalid LP values are
excluded. LD mode accepts biallelic sequence variants on chromosomes 1-22 and X,
including normalized indels. Multiallelic and symbolic records are skipped.

GWAS and reference coordinates must use the same genome assembly. Normalize
variant representations and split multiallelic records before analysis where
needed. Rebuild indexes whenever input files change.

### Alzheimer's disease on chromosome 19

Find Alzheimer’s disease loci on chromosome 19 using the genome-wide GRCh38
summary statistics from Bellenguez et al. (2022) and a European LD reference.
Download the full study and genome sequence once. Only the chromosome 19 LD
panel is needed for this regional analysis.
[GWAS Catalog GCST90027158](https://www.ebi.ac.uk/gwas/studies/GCST90027158) ·
[Study](https://doi.org/10.1038/s41588-022-01024-z)

Run the blocks in order in Bash with `bcftools`, `+gwaslocusfinder`, `+munge`,
`curl`, `gzip`, and `awk`. Downloads total about 2 GiB. Allow 8 GiB of free
disk space for the inputs, sorting, and outputs. Either build of
gwaslocusfinder can use the BCF panel. [+munge installation](https://github.com/freeseek/score#installation)

Download the study from the GWAS Catalog and the GRCh38 no-ALT analysis set
with UCSC chromosome names from NCBI.

```bash
set -euo pipefail
mkdir -p ad-loci-example/tmp
cd ad-loci-example
export TMPDIR="$PWD/tmp"

gwas=GCST90027158_buildGRCh38.tsv.gz
curl -fsSL --retry 3 -o "$gwas" \
  "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/$gwas"

fasta=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
fasta_url=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids
curl -fsSL --retry 3 -o "$fasta.gz" "$fasta_url/$fasta.gz"
gzip -dc "$fasta.gz" > "$fasta"
```

Convert the entire study to GWAS-VCF. `+munge` orients alleles against the
FASTA and creates its `.fai` index if needed. The pipeline retains passing
biallelic records and writes a coordinate-sorted, indexed `ad.bcf` with one
phenotype, `AD_2022`. This BCF can be reused for any chromosome.

```bash
bcftools +munge -c SSF -f "$fasta" -s AD_2022 -Ou "$gwas" |
  bcftools view -f PASS,. -m2 -M2 -Ou |
  bcftools sort -T "$TMPDIR/" -Ob -o ad.bcf -W
```

For this release, `ad.bcf` contains 19,811,681 records, including 435,065 on chromosome 19.

Download the phased HGDP+1kGP v2 reference BCF and its CSI index. Intersect
its sample names with the 1000 Genomes Phase 3 `EUR` labels to create the
reference `--keep` file. The intersection contains **492 individuals**.
[HGDP+1kGP reference publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC11216312/)

```bash
panel=hgdp1kgp_chr19.filtered.SNV_INDEL.phased.shapeit5.bcf
panel_url=https://storage.googleapis.com/gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2
curl -fsSL --retry 3 -o "$panel" "$panel_url/$panel"
curl -fsSL --retry 3 -o "$panel.csi" "$panel_url/$panel.csi"
curl -fsSL --retry 3 -o phase3.panel \
  https://1000genomes.s3.amazonaws.com/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

bcftools query -l "$panel" > ref.samples
awk 'NR == FNR {if ($3 == "EUR") eur[$1] = 1; next} $1 in eur {print $1}' \
  phase3.panel ref.samples > eur.keep
wc -l eur.keep
```

Use the bcftools region option `-r chr19` to select chromosome 19 from the
indexed, genome-wide `ad.bcf`. Run LD clumping and aggregation with the default
significance and LD thresholds and a minimum reference MAF of 1% in the 492
individuals. This panel supplies phased `GT`. HTSlib 1.25 can print a header
warning for its phase-confidence `PP` field while reading it. The plugin uses `GT` here.

```bash
bcftools +gwaslocusfinder -r chr19 -Ob -W -o ad.index.bcf ad.bcf -- \
  --vcf-ref "$panel" --keep eur.keep --maf 0.01 \
  --leads ad.leads.bcf --independent ad.independent.bcf --reports-prefix ad

bcftools query -f '%CHROM\t%POS\t%ID[\t%LF_LOCUS\t%LF_START\t%LF_END]\n' ad.index.bcf
```

The query prints the four index variants, their locus numbers, and their
inclusive bounds on GRCh38. The accompanying outputs contain 65 leads and
184 independent significant variants. For example, `ad.loci.tsv` assigns
58 leads to locus 3, whose index is `rs117310449`.

```text
chr19	1050875	rs12151021	1	1009366	1093669
chr19	1854152	rs12972720	2	1808570	1891993
chr19	44890259	rs117310449	3	43752735	45597073
chr19	54265512	rs12984029	4	54251846	54283071
```

**Using PGEN.** With PLINK 2 and a PGEN-enabled plugin, convert the same
reference and keep the selected individuals at import. This produces the same
index variants, leads, independent variants, and five result tables.

```bash
plink2 --bcf "$panel" --keep eur.keep --make-pgen --out chr19.eur
bcftools +gwaslocusfinder -r chr19 -Ob -W -o ad.pgen.index.bcf ad.bcf -- \
  --pfile chr19.eur --maf 0.01 \
  --leads ad.pgen.leads.bcf --independent ad.pgen.independent.bcf --reports-prefix ad.pgen
```

To annotate the index positions in `ad.index.bcf` with nearby genes, use a
GRCh38 GTF or GFF3 with the optional bcftools plugin [+genecontext](https://github.com/bgorm/genecontext).

### Run an analysis

For a reference fileset `reference.pgen`, `reference.pvar`, and
`reference.psam`:

```sh
bcftools +gwaslocusfinder -Ob -o trait.index.bcf -W trait.bcf -- \
  --pfile reference
```

This analyzes all phenotype columns and writes the union of their index
variants to `trait.index.bcf`, with a CSI index. The default settings select
signals at P <= 5e-8, clump at r2 thresholds of 0.6 and 0.1, evaluate LD within
1 Mb, and merge locus intervals separated by at most 250 kb. See
[Locus definition](#locus-definition) for how these settings interact.

General options and the GWAS input go **before `--`**. Reference-panel,
clumping, and optional-output options go **after `--`**. For example:

```sh
bcftools +gwaslocusfinder -s TRAIT -Ob -o trait.index.bcf -W trait.bcf -- \
  --pfile reference --p1 5e-8 --p2 0.05 \
  --indep-r2 0.6 --lead-r2 0.1 --ld-kb 1000 --merge-kb 250 \
  --reports-prefix trait.loci
```

Inspect index positions and locus boundaries with bcftools:

```sh
bcftools query -f '%CHROM\t%POS\t%ID[\t%SAMPLE:%LF_LOCUS:%LF_START-%LF_END]\n' \
  trait.index.bcf
```

### Select phenotypes and filter GWAS records

`-s/--samples` and `-S/--samples-file` select **GWAS phenotype columns**;
`--keep` selects **reference individuals**. All GWAS phenotypes are analyzed
by default. A sample file contains one phenotype name per line; prefix a list
or file with `^` to exclude names. Unknown names are errors unless
`--force-samples` is set.

To analyze two phenotypes using only GWAS SNVs with FILTER `PASS` or `.`, and
restrict LD estimation to a reference subset:

```sh
bcftools +gwaslocusfinder -s BMI,HEIGHT -f PASS,. -i 'TYPE="snp"' \
  -Oz -o traits.index.vcf.gz -W traits.bcf -- \
  --pfile reference --keep eur_unrelated.keep
```

`-i` and `-e` use standard bcftools expressions to select or reject whole GWAS
records, as in `bcftools view`. They are mutually exclusive. Expressions see
the original sample columns before `-s/-S` subsetting. `-f` optionally selects
FILTER labels; no FILTER selection is implicit.

`--gwas-maf` filters significant signals separately for each phenotype using
its `FORMAT/AF` (`MAF = min(AF, 1-AF)`). A variant can pass for one phenotype
and fail for another. This filter is disabled by default.

Rejected GWAS records cannot return as reference-only candidates. Variants
absent from the GWAS can still extend locus boundaries; add `--gwas-only` to
exclude them. To restrict the reference itself to SNVs while keeping its
non-GWAS variants eligible, prepare a reference containing only SNVs.

`-r/-R` and `-t/-T` restrict significant-signal discovery. Candidate expansion
can extend beyond those regions. Region strings use 1-based inclusive
coordinates; BED files use 0-based half-open coordinates. The default overlap
modes are record overlap for regions and POS overlap for targets.

## Reference panels

### PGEN files and shards

PGEN references require the default build with `WITH_PGEN=1`.

`--pfile` accepts either a fileset prefix or a `.pgen` filename. Each PGEN
requires its matching `.psam` and exactly one of `.pvar`, `.pvar.zst`, or
`.pvar.gz`. PGEN genotypes and retained dosages are read directly with pgenlib.

For a sharded panel, create a list of prefixes or PGEN filenames. Paths are
relative to the list's directory:

```text
reference/chr1.pgen
reference/chr2.pgen
reference/chr3.pgen
```

```sh
bcftools +gwaslocusfinder -Ob -o trait.index.bcf -W trait.bcf -- \
  --pfile-list references.list --keep eur_unrelated.keep \
  --maf 0.01 --min-call-rate 0.95
```

`--pfile` is also repeatable. All shards must have the same FID/IID pairs in
the same order. Shards may split a chromosome; LD is evaluated across shard
boundaries. Panels split by samples are not supported. Duplicate variants
within queried windows generate warnings; the first occurrence in panel-list
and variant-row order is retained, without comparing or merging genotypes.

`--keep` accepts IID or FID IID per line, with an optional header beginning
with `#`. Unknown or ambiguous IDs are errors. All selected individuals,
including nonfounders, contribute to LD; supply the intended ancestry and
unrelated-sample subset explicitly.

PGEN dosages take precedence over hardcalls. To preserve dosages when importing
an autosomal VCF with PLINK, request the dosage field explicitly:

```sh
plink2 --vcf reference.vcf.gz dosage=DS --make-pgen --out reference
```

Use `dosage=HDS` for Minimac-style phased dosage input. PLINK does not import
dosages by default; gwaslocusfinder can use only those retained in the PGEN.
See the [PLINK import documentation](https://www.cog-genomics.org/plink/2.0/input)
for chromosome X import requirements.

### VCF/BCF references

Use `--vcf-ref` for an indexed reference, or `--vcf-ref-list` for shards:

```sh
bcftools +gwaslocusfinder -Ob -o trait.index.bcf -W trait.bcf -- \
  --vcf-ref reference.bcf --keep eur_unrelated.keep

bcftools +gwaslocusfinder -Ob -o trait.index.bcf -W trait.bcf -- \
  --vcf-ref-list references.list --keep eur_unrelated.keep
```

VCF/BCF references have the same local-file and index requirements as the GWAS.
`--vcf-ref` is repeatable, and list paths are relative to the list's directory.
VCF and BCF shards can be mixed, but PGEN and VCF/BCF cannot be mixed in one run.
All VCF/BCF shards must have identical sample names in identical order.
For `--keep`, use VCF sample names as IIDs; use FID `0` in a two-column file.

For each sample, the reader uses **DS, then the sum of two HDS values, then GT**.
Missing DS falls back to HDS; a missing or incomplete HDS pair falls back to GT.
GT-only, DS-only, and HDS-only records are supported. `--hardcalls` ignores
DS/HDS and requires GT in each accessed record. GP is not decoded.

### Chromosome X and diploid reference encoding

Reference values must use a diploid encoding on **all supported chromosomes,
including X**. `X`, `chrX`, and `23` are recognized as chromosome X.

| Reference field | Required encoding |
| --- | --- |
| GT | Two alleles when called; encode haploid REF/ALT calls as 0/0 or 1/1 |
| DS | One Float ALT dosage per sample, on the 0-2 scale |
| HDS | Two Float haplotype ALT dosages, each in [0,1] |
| PGEN | Stored hardcalls and total dosages already use 0-2 encoding |

DS is used unchanged. Already diploid-scaled imputation dosages need no
conversion; haploid 0-1 dosages must be converted before use. The range of
observed values alone cannot identify the intended scale. Called GT and
nonmissing HDS are validated even when DS is available; a single nonmissing
haploid HDS value is rejected.

For a mixed-ploidy VCF/BCF reference, `+fixploidy` can duplicate haploid GT alleles:

```sh
bcftools +fixploidy -Ob -o reference.diploid.bcf reference.bcf -- -f 2
bcftools index -f reference.diploid.bcf
bcftools +gwaslocusfinder -Ob -o trait.index.bcf -W trait.bcf -- \
  --vcf-ref reference.diploid.bcf --hardcalls
```

`+fixploidy` changes GT only; it does not rescale DS or expand HDS. The example
uses hardcalls. Omit `--hardcalls` when the dosage fields also meet the encoding
requirements. PGEN needs no additional scaling.

LD gives equal weight to each selected sample. MAF uses mean dosage divided
by two, without sex-specific chromosome counting. PSAM `SEX` is not used.
Pseudoautosomal variants labeled X follow the same model; `PAR1`, `PAR2`, and
`XY` labels are not automatically remapped to X.

## Locus definition

The independent-variant, lead-variant, and genomic-locus hierarchy follows the
[FUMA locus-discovery model](https://fuma.ctglab.nl/tutorial#riskloci).
Gwaslocusfinder applies the following steps separately to each phenotype:

1. Select variants with P <= `--p1` (default 5e-8), match them to the reference,
   and remove those failing reference QC. Process the remainder from strongest
   to weakest association. Ties use input contig order, position, and allele key.
2. Retain **independent significant variants** by greedy clumping at
   `--indep-r2` (0.6). A variant is retained when its r2 with every stronger
   retained variant within the LD window is below this threshold.
3. Clump the independent variants again at `--lead-r2` (0.1) to select **leads**.
4. Expand a block around each independent variant using neighbors with
   r2 >= `--indep-r2`. GWAS-tagged neighbors require P < `--p2` (0.05);
   reference-only neighbors are eligible unless `--gwas-only` is set. Each
   block includes the independent variant itself.
5. Join independent-variant blocks connected at r2 >= `--lead-r2`, including
   transitive chains. Merge the resulting intervals when their gap is at most
   `--merge-kb` (250). The strongest variant in each locus is its **index**.

In LD mode, **index variants are a subset of leads, which are a subset of
independent significant variants**. Independence refers to clumping within the
LD window; this procedure does not perform conditional association testing.

`--ld-kb` (1000) limits every pairwise LD comparison. More distant pairs are
not evaluated. Block boundaries are the extreme candidate positions, with no
padding by the search radius. Reference-only variants can therefore extend a
block and cause neighboring blocks to merge, but they cannot become index or
lead variants because they have no GWAS association statistic.

`--merge-kb 0` prevents merging across positive gaps. Overlapping blocks and
blocks connected by LD can still form one locus, so index output can contain
fewer records than lead output. Merging changes lead locus assignments, not
which lead variants are retained. For a single effective clumping threshold,
set `--indep-r2` and `--lead-r2` to the same value. `--p2` must be at least
`--p1`; it controls candidate blocks, not a second significance threshold.

### Variant matching and LD

Variants are matched by chromosome, position, and alleles, not by rsID.
A leading `chr` is ignored. Single-base allele pairs may be reversed because
squared correlation is unchanged; REF/ALT order is preserved for indels and
other multi-base alleles. The tool does not complement strands or normalize
variants. Canonical SNP keys do not specify the effect allele.

LD is squared Pearson correlation of unphased ALT dosages using pairwise-complete
samples. Phase is ignored. `--hardcalls` uses hardcall allele counts instead.
Reference `--maf` and `--min-call-rate` are calculated in the selected sample
subset; monomorphic or constant-dosage variants are always excluded.
The initial `nearby variants` count is before these checks. The final
`Reference QC` line reports distinct variants evaluated, passed, filtered out,
and not evaluated (genotypes not needed) across all selected phenotypes.
These counts also appear as `reference_variants_qc_evaluated`,
`reference_variants_qc_passed`, `reference_variants_qc_filtered`, and
`reference_variants_qc_not_evaluated` in `.run.tsv`. Passed + filtered out +
not evaluated equals `reference_variants_in_windows`.

QC runs only when a variant's genotypes are needed. Variants outside the
analyzed candidate blocks, or skipped by candidate eligibility rules, may
remain unevaluated. The passing count therefore describes evaluated variants;
it is not an estimate for the unevaluated part of the panel. Counting adds no
extra genotype scan.

Missing reference signals are skipped and reported by default; use
`--missing-ref error` to require their presence.

At least `--min-pair-samples` complete pairs (default 3) are required. Undefined
LD between potential significant/independent signals is an error; undefined
candidate pairs are skipped. Stored PGEN dosage precision is 1/16384, whereas
VCF/BCF dosages retain Float precision, so conversion can affect decisions
very close to an r2 cutoff.

The implementation uses the FUMA-style hierarchy with the supplied reference
and these explicit LD rules. Agreement with a particular FUMA analysis depends
on its reference, sample selection, filters, and LD calculation. Functional
annotation and gene mapping are outside the scope of this tool.

### Distance mode

Without a reference panel, use `--distance-kb`:

```sh
bcftools +gwaslocusfinder -Ob -o trait.index.bcf -W trait.bcf -- \
  --distance-kb 1000
```

Consecutive significant positions separated by at most this distance are
joined transitively. The strongest variant in each interval is its index.
This mode makes no LD-independence claim; `--leads` and `--independent` require
an LD reference panel.

## Outputs

### Index, lead, and independent variants

The main output is the **set of index variants**, one per locus and phenotype.
Use `-o` to name the file; the plugin defaults to plain VCF on standard output.
Request additional selections from the same analysis with `--leads` and
`--independent`:

```sh
bcftools +gwaslocusfinder -Oz -o trait.index.vcf.gz -W=tbi trait.bcf -- \
  --pfile reference \
  --leads trait.leads.vcf.gz --independent trait.independent.vcf.gz \
  --reports-prefix trait.loci
```

`-O` controls every variant output: `v` for plain VCF, `z` for compressed VCF,
`b` for compressed BCF, or `u` for uncompressed BCF. A digit sets the compression
level, for example `-Oz9`. File extensions do not select the format.
`-W` creates CSI indexes for named compressed outputs; `-W=tbi` requests TBI
for compressed VCF. Standard output is not indexed. At most one variant stream
can use stdout; name the main file when an optional output uses `-`.

For LD-clumped leads without reducing each locus to one index variant:

```sh
bcftools +gwaslocusfinder -Ob -o /dev/null trait.bcf -- \
  --pfile reference --indep-r2 0.1 --lead-r2 0.1 --leads trait.leads.bcf
bcftools index trait.leads.bcf
```

The main stream goes to `/dev/null`. Omit `-W` in this case because it applies
to every named variant output, including the device; index the leads separately.
Each named output, report, and index sidecar needs its own path. The plugin
checks these paths against the GWAS, the reference, sample and interval files,
reference lists, and file operands in bcftools filters before writing.
Existing output files are overwritten directly.

### Phenotypes and locus annotations

Each output contains the union of variants selected by the requested phenotypes.
A record is written once even if several phenotypes select it. New annotations
are stored in FORMAT:

| Field | Meaning |
| --- | --- |
| LF_LOCUS | Locus number within the phenotype |
| LF_START | Parent locus start, 1-based inclusive |
| LF_END | Parent locus end, 1-based inclusive |
| LF_INDEX_POS | Position of the strongest variant in the parent locus, 1-based |

All four FORMAT fields have `Number=1,Type=Integer`. An input from an earlier
run can supply those definitions. The plugin writes fresh values for this run
and reports an error if a reserved FORMAT field has an incompatible definition.

Locus IDs start at 1 separately for each phenotype and are shared across the
index, lead, and independent outputs from a run. Join by phenotype and
`LF_LOCUS`, or by `(SAMPLE, LOCUS)` in the TSVs. Bounds and index positions are
on the record's chromosome.

Where a phenotype did not select a variant for a particular output, **all of
its FORMAT values are missing**. Selected combinations retain their original
GWAS values. For example:

| Variant | A: LP | A: LF_LOCUS | B: LP | B: LF_LOCUS |
| --- | --- | --- | --- | --- |
| rs1 | 30 | 1 | . | . |
| rs2 | . | . | 40 | 1 |

All existing INFO fields are preserved. Records retain their original positions,
alleles, and IDs and are coordinate sorted. Locus annotations do not change a
variant's span for ordinary VCF region queries.

### Tabular reports

`--reports-prefix PREFIX` adds six tab-separated reports:

| Suffix | Contents |
| --- | --- |
| .loci.tsv | Locus boundaries, index variant, and independent/lead/candidate counts |
| .independent.tsv | Independent significant variants, assigned leads, and parent loci |
| .leads.tsv | Lead variants and parent loci |
| .candidates.tsv | One row per independent-variant/candidate LD relationship |
| .excluded.tsv | Significant variants excluded for missing reference, QC, or unsupported chromosome |
| .run.tsv | Parameters, input paths, selected samples, versions, and run counts |

Result tables include `SAMPLE`; locus tables use `LOCUS`. The lead and independent
tables include `LOCUS_START`, `LOCUS_END`, `INDEX_KEY`, `INDEX_ID`, and `INDEX_POS`.
Candidate REF/ALT values come from the reference. The run report is written after
all requested outputs succeed. Empty analyses produce valid empty variant files
and header-only result tables. In distance mode, LD tables are header-only and
LD counts are missing.

## Clumping and thresholding for polygenic scores

Combine gwaslocusfinder with [`+munge` and `+score`](https://github.com/freeseek/score)
to convert summary statistics, retain LD-clumped weights, and compute polygenic
scores across a P-value grid. The weights are the original marginal effects in
`FORMAT/ES`; a separate weight-estimation step is not needed for this C+T method.
Target genotypes are required in addition to the GWAS and LD reference panel.

### Prepare and harmonize the inputs

The example uses a single phenotype and autosomal SNVs on GRCh38. Prepare:

- `sumstats.tsv.gz`: quality-controlled summary statistics with positions,
  effect/other alleles, signed effects, and P-values, in GWAS-SSF format.
- `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna` and its `.fai`:
  the GRCh38 reference sequence.
- `reference.pgen/.pvar/.psam`: an LD panel representing the intended ancestry.
- `target.qc.bcf` and its index: sorted, biallelic target genotypes with ALT
  dosages in `FORMAT/DS`.

All inputs must use the same assembly and compatible allele representations.
The GWAS and target must also have identical contig names for exact matching.
Resolve duplicate records and strand ambiguity before clumping. FASTA matching
alone does not resolve ambiguous A/T or C/G strand conventions. For different
assemblies, use [`+liftover`](https://github.com/freeseek/score#liftover-vcfs),
then normalize, sort, and index before intersecting the variant sets.

Convert summary statistics to GWAS-VCF:

```sh
set -euo pipefail
bcftools +munge -c SSF \
  -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  -s TRAIT -Ou sumstats.tsv.gz | \
  bcftools view -f PASS,. -m2 -M2 -Ou | \
  bcftools sort -Ob -o trait.bcf -W
```

For other column layouts, use the appropriate preset or `-C colheaders.tsv`.
Map the effect allele to `A1`, beta/log-odds coefficients to `BETA`, or raw odds
ratios to `OR`. With `-f`, `+munge` orients ES to ALT, changing its sign when
needed, and converts P-values to LP. The pipeline excludes unresolved or
mismatching reference-allele calls. Check the conversion log and ES/LP values;
gwaslocusfinder cannot supply missing effect estimates.

Restrict the GWAS to alleles available in the target **before clumping**:

```sh
bcftools isec -c none -n=2 -w1 -Ob -o trait.shared.bcf -W \
  trait.bcf target.qc.bcf
```

This matches chromosome, position, REF, and ALT and retains the GWAS records.
A stronger variant absent from the target can otherwise suppress an available
proxy during clumping. For multiple target cohorts, define a consistent variant
pool. Missing LD-reference matches are skipped by default; use
`--reports-prefix` to inspect exclusions or `--missing-ref error` to reject them.

### Clump and score

Use **leads** as the scoring weights. The main index output can omit additional
LD-clumped leads in a merged locus. Here both clumping thresholds are 0.1, with
LD evaluated within 250 kb:

```sh
bcftools +gwaslocusfinder -s TRAIT \
  -i 'TYPE="snp" && FMT/ES!="."' \
  -Ob -o /dev/null trait.shared.bcf -- \
  --pfile reference --p1 1 --p2 1 \
  --indep-r2 0.1 --lead-r2 0.1 --ld-kb 250 \
  --gwas-only --leads trait.clumped.bcf
bcftools index -f trait.clumped.bcf

bcftools +score --use DS --counts \
  --q-score-thr 5e-8,1e-6,1e-4,0.001,0.01,0.05,0.1,0.5,1 \
  -o trait.prs.tsv target.qc.bcf trait.clumped.bcf
```

`--p1 1` retains variants for thresholds up to 1; `--p2 1` satisfies the
requirement that p2 >= p1. The input must include variants up to the largest
intended cutoff. Clump once at that cutoff for a fixed LD threshold, radius,
and eligible variant set: weaker variants cannot displace stronger retained
variants. Changing the LD settings requires another clumping run.

The expression excludes missing effects before clumping. `--gwas-only` limits
candidate expansion, while locus merging changes annotations without changing
the retained leads. For sharded or VCF/BCF LD panels, substitute
`--pfile-list references.list` or `--vcf-ref reference.bcf`. The **target** for
`+score` must be VCF/BCF; export a PGEN target with its dosages when necessary.

For each individual and threshold, `+score` sums `ES * ALT dosage` over retained
variants passing the cutoff with usable genotypes. It uses ES and LP, not the
`LF_*` annotations. The output contains scores and contributing-variant counts
for each threshold. `-s/-S` select GWAS phenotypes in gwaslocusfinder but target
individuals in `+score`. Exact threshold membership can depend on stored LP
precision and the scorer's cutoff comparison.

With `+score --use DS`, missing DS values are skipped without mean imputation
or fallback to GT. Use `--use GT` for hardcall scores, or prepare missing dosages
according to the intended scoring policy. Inspect matched-marker logs and
`--counts` for lost variants and uneven coverage.

Choose P-value and LD settings in tuning data, then evaluate the fixed choice
in held-out data. Ancestry matching, discovery-sample overlap, covariate
adjustment, score standardization, and calibration need to be addressed in the
analysis design. These commands calculate raw scores; they do not select a
model or evaluate predictive performance.

## Performance and limitations

`--threads` controls HTSlib I/O; LD calculations are currently serial.
`--cache-mb` bounds the decoded genotype cache, not total process memory.
Reference metadata is sorted into temporary binary files, with a sparse
coordinate index and a small block cache in RAM. Phenotype P-values and QC
annotations are also stored on disk. Candidate relationships are spooled only
when TSV reports are requested. This requires no additional library.

Use `--temp-dir` to choose scratch storage with sufficient free space:

```sh
bcftools +gwaslocusfinder -Ob -o index.bcf trait.bcf -- \
  --pfile reference --temp-dir /path/to/scratch
```

The default is `$TMPDIR` when set, otherwise the current directory. Choose a
disk-backed directory: a RAM-backed `/tmp` still consumes memory. A NAS share
is supported, although its I/O speed affects runtime. Temporary files are
removed on completion and on reported errors; an interrupted process may leave
its `gwaslocusfinder.*` directory for manual removal. Sorting can temporarily
require space for both input runs and merged files; disk use increases with
reference coverage, phenotype count, and optional candidate reports.

PVAR is scanned once to preserve original PGEN row numbers. PGEN's own
whole-file indexes, the sparse coordinate index, and retained significant GWAS
records use RAM across LD windows. Full reference metadata and candidate
relationships are stored in temporary files.

The second GWAS pass uses the index to visit merged LD windows containing
reference variants. It can read beyond discovery regions to tag candidates at
both ends of each window. The run report records the GWAS rows visited
(`tagging_records`), time spent tagging (`tagging_seconds`), and time spent
expanding candidate blocks (`candidate_expansion_seconds`).

VCF/BCF reference scans use indexed signal windows. The metadata pass skips
FORMAT parsing; genotype validation occurs when records are fetched for LD.
Nearby forward genotype reads reuse the reader within the same signal window.
Text VCF still requires reading and decompressing whole selected records,
including sample columns, so broad window coverage can involve substantial
I/O. BCF avoids text parsing and is generally preferable for repeated analyses.

Broad P-value cutoffs can substantially increase work; candidate blocks and
aggregated loci are still constructed when only lead variants are needed for
scoring.

LD mode excludes chromosomes other than 1-22 and X, reporting significant
records as `unsupported_chromosome`. Sex-specific weighting, automatic PAR
remapping, and multiallelic allele-dosage decomposition are not implemented.
For measured wall time, memory, scratch space, and genotype reads on real
GWAS inputs, see the [performance benchmarks](benchmarks/README.md).
Agreement with a particular FUMA analysis requires separate validation.

## License

The gwaslocusfinder source code is distributed under the [MIT License](LICENSE).
HTSlib and optional pgenlib retain their own licenses; see
[THIRD_PARTY.md](THIRD_PARTY.md) for sources and dependency notices.
