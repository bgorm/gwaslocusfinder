# Reading a sample subset from PGEN

The baseline asks pgenlib for every sample, then walks the full row to collect
the people listed in `--keep`. The candidate passes their sample mask to
pgenlib once per shard and converts the returned genotypes and dosages to
floats in panel order. When all samples are selected it uses pgenlib's
full-panel read and converts its output directly. The measured commits were
`50507a2` and `0972478`, respectively.
[Baseline source](https://github.com/bgorm/gwaslocusfinder/commit/50507a2adec2a68d3ca218bfc1557705327b8414) ·
[Candidate reader](../../src/reference.hpp)

With **396 of 3,163 people** selected, preparing genotype rows for the EA4
chromosome 22 analysis took **1.240 → 0.310 s**. The complete analysis took
**6.012 → 5.014 s**, a 16.6% reduction. All six paired runs finished sooner.
On a synthetic panel with stored fractional dosages, time to prepare rows
fell from **0.600 to 0.200 s** and total time fell by **44.7%**.

| Workload and samples kept | Genotype reads | pgenlib (s) | Ready row (s) | Wall (s) | Individual runs |
| --- | ---: | ---: | ---: | ---: | --- |
| EA4, 1 Mb chr22, 396 | 60,411 | 0.013 → 0.014 | 0.227 → 0.070 | 1.674 → 1.465 | [JSON](measurements/regional-396.json) |
| EA4, chr22, 396 | 361,202 | 0.085 → 0.093 | 1.240 → 0.310 | 6.012 → 5.014 | [JSON](measurements/chr22-396.json) |
| Synthetic, dosages, 396 | 30,000 | 0.044 → 0.114 | 0.600 → 0.200 | 0.873 → 0.483 | [JSON](measurements/dosage-396.json) |
| Same synthetic panel, `--hardcalls`, 396 | 30,000 | 0.029 → 0.042 | 0.367 → 0.101 | 0.722 → 0.451 | [JSON](measurements/hardcalls-396.json) |
| EA4, 1 Mb chr22, 3,005 | 63,037 | 0.017 → 0.020 | 0.368 → 0.239 | 2.166 → 1.995 | [JSON](measurements/regional-3005.json) |
| EA4, 1 Mb chr22, all 3,163 | 63,053 | 0.017 → 0.017 | 0.385 → 0.247 | 2.120 → 2.070 | [JSON](measurements/regional-3163.json) |
| Synthetic, dosages, 3,005 | 30,000 | 0.044 → 0.110 | 0.620 → 0.628 | 1.681 → 1.615 | [JSON](measurements/dosage-3005.json) |

Each cell shows the baseline followed by the candidate. Applying a sample
mask inside pgenlib took extra time on the dosage panel. With 396 samples,
it reduced the packed dosages passed back to the plugin from 34.57 million
to 4.33 million. Converting just the selected samples accounted for the gain.
With 3,005 people, ready-row time was essentially even (three paired wins
each). In every workload, both readers returned byte-identical genotypes in
the same order and produced byte-identical variant and table outputs.

## What the timers measure

- **pgenlib** measures time inside `PgrGetD` (or `PgrGet` for `--hardcalls`),
  summed across genotype reads.
- **Ready row** includes that call, allocation of the plugin's float vector,
  dosage or hardcall conversion, and sample selection. This is the time to
  deliver the genotype values used by QC and LD.
- **Wall** measures the complete uninstrumented bcftools command, including
  discovery, QC, LD, locus aggregation, three BCF outputs and the reports.

Each workload used six paired runs for each timer configuration, alternating
which reader ran first. The tables report medians across the six runs.
Verification runs warmed the filesystem cache before timing. The dedicated
profile builds add three steady-clock reads per genotype read. Both readers
used the same host, compiler flags, inputs and analysis parameters. Measurements
were sequential on Linux/WSL2 x86-64, AMD Ryzen 9 5900XT, with inputs, output
and scratch on NFS, GCC 15.2 `-O3 -g`, bcftools/HTSlib 1.25, zero HTSlib
worker threads, `--maf 0.05`, and `--cache-mb 128`.
[Build and binary checksums](measurements/builds.json)

The EA4 analysis used the 1,093,149-variant chromosome 22 PGEN shard from
the regional benchmark. The rows read contained hardcalls only. The synthetic
panel supplied fractional dosages, missing values and a single significant
variant whose window covers all 30,000 rows. The 396-person selection kept
every eighth of 3,163 samples in panel order. The 3,005-person selection
dropped every twentieth sample.
[EA4 input preparation](../README.md#reproduce-the-measurements)

In a separate verification run per reader and workload, the probe hashed
the stream of PGEN row numbers and every returned float in sample order,
including missing values. Each timed run also checked the number of genotype
reads, reference samples, variants passing QC, loci, and SHA-256 of all three
BCFs and five result tables. The JSON files include the individual runs and
digests. The dosage counters describe the packed arrays pgenlib returned, so
they count the full panel in the baseline and only selected samples in the
candidate. [Timing and verification code](measure.py) · [Probe](pgen_probe.hpp)

## Reproduce the comparison

From a clone containing this benchmark, build both pinned revisions against
one set of matching headers and pgenlib sources. These commands require
bcftools 1.25 on the path. Use a new directory with room for the reference
panel and scratch.

```bash
set -euo pipefail
method="$PWD/benchmarks/pgen-subset"
study="$PWD/../pgen-subset-benchmark"
mkdir "$study"
git worktree add --detach "$study/baseline" 50507a2adec2a68d3ca218bfc1557705327b8414
git worktree add --detach "$study/candidate" 0972478d05a57f1830b8b488374a35c2831c3a7c
(cd "$study/candidate" && bash scripts/bootstrap.sh)
cd "$study"
export PLUGIN_HTSLIB="${PLUGIN_HTSLIB:-$study/candidate/vendor/htslib-1.25}"
export PGENLIB="${PGENLIB:-$study/candidate/vendor/plink-ng-4544f92e374e784a6ffd97fa88239962b81f2b62/2.0/include}"
make -C baseline -j4 plugin
for reader in baseline candidate; do
  python3 "$method/instrument.py" "$reader" "profile-$reader"
  make -C "profile-$reader" -j4 plugin
done
```

Prepare the EA4 GWAS, chromosome 22 shard and every-eighth keep list as in
the regional benchmark linked above. To generate the synthetic panel, the
measured runs used PLINK v2.0.0-b.1.3-devLM (26 Oct 2026). Its dummy marker
positions start at zero. Shift them to one-based coordinates and write a
GWAS record matching the middle marker.

```bash
mkdir synthetic
plink2 --dummy 3163 30000 0.03 acgt dosage-freq=0.5 phase-freq=0.5 \
  --seed 32683 --threads 1 --make-pgen --out synthetic/panel
awk 'BEGIN {OFS="\t"} !/^#/ {$2+=1} {print}' synthetic/panel.pvar > synthetic/shifted.pvar
mv synthetic/shifted.pvar synthetic/panel.pvar
{
  printf '%s\n' '##fileformat=VCFv4.2' '##contig=<ID=1,length=1000000>' \
    '##FORMAT=<ID=LP,Number=A,Type=Float,Description="Negative log10 P">'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttrait\n'
  awk 'BEGIN {OFS="\t"} $3=="snp15000" {print $1,$2,$3,$4,$5,".","PASS",".","LP",30}' synthetic/panel.pvar
} | bcftools view --no-version -Ob -W -o synthetic/gwas.bcf
```

Set the three EA4 input paths and choose a new results directory. The driver
creates the denser keep files, runs the seven workloads, and writes each
completed pair to JSON. It streams the verification data to SHA-256.

```bash
python3 "$method/measure.py" "$study" /path/to/ea.b38.vcf.gz \
  /path/to/runs/panel22 /path/to/runs/every8.keep \
  synthetic/panel synthetic/gwas.bcf measurements
```
