# Performance on GWAS-VCF inputs

Whole-genome locus discovery on the SSGAC EA4 GWAS took **18 min 34 s**, peaked
at **1.81 GiB RSS**, and used **4.09 GiB** of scratch at the largest sampled
point. The GRCh38 input contained 10,983,047 records. The PGEN reference had
75,193,455 rows and 3,202 individuals. With `--maf 0.05`, the plugin found
619 loci, 1,176 leads, and 4,437 independent significant variants.

| Workload | GWAS rows visited in tagging pass | Reference variants in LD windows | Genotype reads | Wall (s) | Peak RSS (MiB) | Sampled peak scratch (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| EA4, whole genome | 4,669,010 | 29,367,313 | 39,121,634 | 1,114.10 | 1,848.53 | 4,192.34 |
| EA4, chromosome 22 | 55,850 | 397,013 | 392,898 | 8.65 | 157.51 | 48.59 |
| EA4, 1 Mb of chromosome 22, v0.6.0 | 10,983,047 | 67,424 | 63,053 | 14.26 | 156.94 | 6.55 |
| EA4, same 1 Mb, indexed tagging | 9,385 | 67,424 | 63,053 | 2.12 | 157.04 | 6.55 |
| EA4, same 1 Mb, 396 reference samples | 9,385 | 67,424 | 60,411 | 1.66 | 127.45 | 6.55 |
| Chromosome 22 BCF, EA4 selected | 141,015 | 397,013 | 315,777 | 7.87 | 155.88 | 47.55 |
| Same BCF, EA4 and GIANT BMI selected | 248,015 | 684,848 | 528,865 | 12.87 | 155.78 | 84.41 |
| EA4, chromosome 22 C+T with `--p1 1` | 147,598 | 1,040,019 | 4,042,420 | 288.24 | 269.06 | 126.71 |

The regional comparison used the full indexed EA4 VCF and the same PGEN shard
for both versions. Each value in those two rows is the median of three runs
per version, alternating versions. The index, lead, and independent BCFs and
all five result tables were byte-identical across the six runs. As a control,
v0.6.0 took 2.15 s and produced the same results when the GWAS file itself was
cropped to the discovery interval plus 1 Mb on both sides.

Other rows show one run each. The chromosome 22 shard contained 1,093,149
variants and the 3,163 individuals shared between the phased HGDP+1kGP v2
panel and the whole-genome reference. The `--keep` row used every eighth
individual in shard sample order. This timing includes QC and LD on the
selected individuals. The PGEN reader still decodes full rows. Genotype reads
count cache misses, including repeated reads of a variant.

Both multi-phenotype cases used one chromosome 22 BCF with 369,218 records.
Records present for either phenotype can change candidate eligibility, so the
one-versus-two comparison uses `-s ea` on that same BCF. In the chromosome-wide
C+T case, candidate expansion accounted for 116.70 of 288.24 seconds and the
lead output contained 3,901 variants. Candidate time includes genotype access,
QC, pairwise LD, and writing qualifying relationships to scratch.

## Reproduce a run

Measurements used a 0.6.1 development build, bcftools and HTSlib 1.25, GCC 15.2
with `-O3`, and Linux/WSL2 on an AMD Ryzen 9 5900XT. Inputs, outputs, and
scratch were on NFS. Each run used the default of zero HTSlib worker threads.
The filesystem cache was left enabled. Disk use and elapsed time depend on the
filesystem, data density, sample count, and cache settings.

The [measurement script](../scripts/benchmark.py) accepts a **new output
directory**, a plugin name or `+path/to/plugin.so`, and ordinary plugin
arguments. It requests all three BCFs and the TSV reports and stores them,
the log, and `metrics.json` in that directory. Wall time covers the entire
bcftools process. Peak RSS comes from the OS. Scratch is allocated space
in the plugin's temporary files sampled at roughly 100 ms intervals. Input
files, final outputs, and data preparation are outside the scratch figure.

The following commands assume indexed GRCh38 GWAS-VCFs for EA4 and GIANT BMI,
a whole-genome 1000 Genomes PGEN fileset, and the chromosome 22 phased BCF from
the [HGDP+1kGP reference release](https://pmc.ncbi.nlm.nih.gov/articles/PMC11216312/).
For a tutorial with direct downloads, see the
[Alzheimer's disease example](../README.md#alzheimers-disease-on-chromosome-19).

```bash
ea=ea.b38.vcf.gz
bmi=bmi.b38.vcf.gz
panel=all_hg38
chr22=hgdp1kgp_chr22.filtered.SNV_INDEL.phased.shapeit5.bcf
mkdir -p runs
plink2 --bcf "$chr22" --keep "$panel.psam" --make-pgen --threads 2 --out runs/panel22
bcftools merge -r chr22 -m none --no-version -Ob -W -o runs/two-traits.bcf "$ea" "$bmi"
awk 'NR>1 && (NR-2)%8==0 {print $1}' runs/panel22.psam > runs/every8.keep

python3 scripts/benchmark.py runs/genome +gwaslocusfinder "$ea" -- \
  --pfile "$panel" --maf 0.05 --cache-mb 512

python3 scripts/benchmark.py runs/regional +gwaslocusfinder \
  -r chr22:29000001-30000000 "$ea" -- --pfile runs/panel22 --maf 0.05 --cache-mb 128
```

For chromosome 22 discovery, change the regional argument to `-r chr22`.
For the sample subset, append `--keep runs/every8.keep` to the regional
command. To reproduce the multi-phenotype comparison, use `runs/two-traits.bcf`
as the GWAS input and compare all samples with `-s ea`. Use a new directory
for each run.

```bash
python3 scripts/benchmark.py runs/ct +gwaslocusfinder \
  -r chr22 -i 'TYPE="snp" && FMT/ES!="."' "$ea" -- \
  --pfile runs/panel22 --maf 0.05 --cache-mb 128 --p1 1 --p2 1 \
  --indep-r2 0.1 --lead-r2 0.1 --ld-kb 250 --gwas-only
```
