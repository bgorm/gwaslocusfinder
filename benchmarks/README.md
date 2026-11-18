# Performance on GWAS-VCF inputs

Whole-genome locus discovery on the SSGAC EA4 GWAS took **17 min 56 s**, peaked at
**1.80 GiB RSS**, and used **4.09 GiB** of scratch at the largest sampled
point. The GRCh38 input contained 10,983,047 records. The PGEN reference had
75,193,455 rows and 3,202 individuals. With `--maf 0.05`, the plugin found
619 loci, 1,176 leads, and 4,437 independent significant variants.

Both PGEN-enabled binaries were built from clean Git checkouts. The baseline
used 0.6.0 at `84cde27` and the indexed runs used 0.6.1 at `1a22388`.
The build manifest records the full Git revisions, plugin binary SHA-256
checksums, and measurement script revision.
[Baseline source](https://github.com/bgorm/gwaslocusfinder/commit/84cde274211342add00216c9523085cb4f60a7d6) ·
[Indexed source](https://github.com/bgorm/gwaslocusfinder/commit/1a22388423ca243adaeb88e03497484c441653f6) ·
[Build manifest](measurements/provenance.json)

| Workload | GWAS rows visited in tagging pass | Reference variants in LD windows | Genotype reads | Wall (s) | Peak RSS (MiB) | Sampled peak scratch (MiB) | Run JSON |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| EA4, whole genome | 4,669,010 | 29,367,313 | 39,121,634 | 1,075.91 | 1,846.92 | 4,188.22 | [JSON](measurements/1a22388/genome/metrics.json) |
| EA4, chromosome 22 | 55,850 | 397,013 | 392,898 | 9.27 | 157.48 | 48.23 | [JSON](measurements/1a22388/chr22/metrics.json) |
| EA4, 1 Mb of chromosome 22, 0.6.0 | 10,983,047 | 67,424 | 63,053 | 14.85 | 156.93 | 6.55 | [1](measurements/84cde27/regional-1/metrics.json) · [2](measurements/84cde27/regional-2/metrics.json) · [3](measurements/84cde27/regional-3/metrics.json) |
| EA4, same 1 Mb, indexed tagging | 9,385 | 67,424 | 63,053 | 2.18 | 156.95 | 6.55 | [1](measurements/1a22388/regional-1/metrics.json) · [2](measurements/1a22388/regional-2/metrics.json) · [3](measurements/1a22388/regional-3/metrics.json) |
| EA4, same 1 Mb, 396 reference samples | 9,385 | 67,424 | 60,411 | 1.69 | 127.89 | 6.55 | [JSON](measurements/1a22388/subset/metrics.json) |
| Chromosome 22 BCF, EA4 selected | 141,015 | 397,013 | 315,777 | 7.78 | 155.33 | 47.41 | [JSON](measurements/1a22388/one-phenotype/metrics.json) |
| Same BCF, EA4 and GIANT BMI selected | 248,015 | 684,848 | 528,865 | 12.97 | 155.62 | 82.94 | [JSON](measurements/1a22388/two-phenotypes/metrics.json) |
| EA4, chromosome 22 C+T with `--p1 1` | 147,598 | 1,040,019 | 4,042,420 | 287.28 | 270.58 | 127.21 | [JSON](measurements/1a22388/ct/metrics.json) |

The regional measurements used the same full GWAS file and PGEN shard.
The two full-panel 1 Mb rows report per-column medians from three runs per
build, alternating versions. The index, lead, and independent BCFs and
all five result tables were byte-identical across the six runs. The baseline
visited the entire GWAS during tagging. Its tagging count comes from the
input index because 0.6.0 did not record this counter. Cropping the GWAS to
the discovery interval plus 1 Mb on each side gave the same outputs and took
2.18 s with the baseline. [Cropped control](measurements/84cde27/cropped-control/metrics.json)

Other rows show one run each. The chromosome 22 PGEN shard contained
1,093,149 variants and the 3,163 individuals shared between the phased
HGDP+1kGP v2 panel and the whole-genome reference. The `--keep` run used
every eighth individual in shard order. Its time includes decoding full PGEN
rows, followed by QC and LD on the selected individuals. Genotype reads count
cache misses, including repeat reads of the same variant.

The one- and two-phenotype runs used the same 369,218-record chromosome 22
BCF. Keeping the input BCF fixed controls which alleles are present in the GWAS.
Selecting both phenotypes expands the LD windows visited during tagging.

For chromosome-wide C+T, candidate expansion took 116.67 of 287.28 seconds
(40.6%) and the lead BCF contained 3,901 variants. Candidate time includes
genotype access, QC, pairwise LD, and scratch writes.

The original indexed measurements came from an uncommitted development build.
All 13 pinned runs reproduced the three BCFs and five result tables from their
corresponding original runs byte-for-byte.
[Original table](https://github.com/bgorm/gwaslocusfinder/blob/1a22388423ca243adaeb88e03497484c441653f6/benchmarks/README.md) ·
[Original individual measurements](measurements/initial)

## Reproduce the measurements

Measurements used bcftools and HTSlib 1.25, GCC 15.2 with `-O3`, and Linux/WSL2
on an AMD Ryzen 9 5900XT. Inputs, outputs, and scratch were on NFS. Each run
used the default of zero HTSlib worker threads. Runs were sequential and the
filesystem cache was left enabled. Wall time covers the complete bcftools
process, including all three variant outputs and the TSV reports. Peak RSS
comes from the OS. Scratch is allocated disk space in the plugin’s temporary
files, sampled about every 100 ms.

From a clone of this repository, build both revisions and use the same
[measurement script](https://github.com/bgorm/gwaslocusfinder/blob/1a22388423ca243adaeb88e03497484c441653f6/scripts/benchmark.py)
for both binaries. Each invocation takes a new output directory.

```bash
set -euo pipefail
baseline_src="$PWD/../gwaslocusfinder-060"
indexed_src="$PWD/../gwaslocusfinder-061"
git worktree add --detach "$baseline_src" 84cde274211342add00216c9523085cb4f60a7d6
git worktree add --detach "$indexed_src" 1a22388423ca243adaeb88e03497484c441653f6
(cd "$baseline_src" && bash scripts/bootstrap.sh)
(cd "$indexed_src" && bash scripts/bootstrap.sh)
baseline="+$baseline_src/build/gwaslocusfinder.so"
indexed="+$indexed_src/build/gwaslocusfinder.so"
measure="$indexed_src/scripts/benchmark.py"
```

Prepare indexed GRCh38 GWAS-VCFs for EA4 (phenotype `ea`) and GIANT BMI,
the whole-genome 1000 Genomes PGEN reference, and the chromosome 22 phased
BCF from the [HGDP+1kGP release](https://pmc.ncbi.nlm.nih.gov/articles/PMC11216312/).
For a smaller example with direct downloads, see the
[Alzheimer’s disease walkthrough](../README.md#alzheimers-disease-on-chromosome-19).

```bash
ea=ea.b38.vcf.gz
bmi=bmi.b38.vcf.gz
panel=all_hg38
chr22=hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.bcf
mkdir -p runs
plink2 --bcf "$chr22" --keep "$panel.psam" --make-pgen --threads 2 --out runs/panel22
bcftools merge -r chr22 -m none --no-version -Ob -W -o runs/two-traits.bcf "$ea" "$bmi"
awk 'NR>1 && (NR-2)%8==0 {print $1}' runs/panel22.psam > runs/every8.keep
bcftools view --no-version -r chr22:28000001-31000000 -Oz -W -o runs/cropped.vcf.gz "$ea"
```

Run the regional comparison, cropped control, subset, and phenotype workloads.

```bash
region=chr22:29000001-30000000
ld=(--pfile runs/panel22 --maf 0.05 --cache-mb 128)
for i in 1 2 3; do
  python3 "$measure" "runs/baseline-$i" "$baseline" -r "$region" "$ea" -- "${ld[@]}"
  python3 "$measure" "runs/indexed-$i" "$indexed" -r "$region" "$ea" -- "${ld[@]}"
done
python3 "$measure" runs/control "$baseline" -r "$region" runs/cropped.vcf.gz -- "${ld[@]}"
python3 "$measure" runs/subset "$indexed" -r "$region" "$ea" -- "${ld[@]}" --keep runs/every8.keep
python3 "$measure" runs/chr22 "$indexed" -r chr22 "$ea" -- "${ld[@]}"
python3 "$measure" runs/one-phenotype "$indexed" -s ea runs/two-traits.bcf -- "${ld[@]}"
python3 "$measure" runs/two-phenotypes "$indexed" runs/two-traits.bcf -- "${ld[@]}"
```

The C+T run clumps chromosome 22 at all P-values. The whole-genome run
uses the default genome-wide significance threshold and a 512 MiB cache.

```bash
python3 "$measure" runs/ct "$indexed" -r chr22 -i 'TYPE="snp" && FMT/ES!="."' "$ea" -- \
  "${ld[@]}" --p1 1 --p2 1 --indep-r2 0.1 --lead-r2 0.1 --ld-kb 250 --gwas-only
python3 "$measure" runs/genome "$indexed" "$ea" -- --pfile "$panel" --maf 0.05 --cache-mb 512
```
