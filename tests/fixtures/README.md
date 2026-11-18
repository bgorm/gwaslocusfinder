# Synthetic LD reference

The seven records and nine samples in `panel.vcf` were constructed for the
public plugin tests. The first record has two alternate alleles. The biallelic
rows include exact quarter dosages, cells with GT and missing DS, one missing
GT, and candidates that are absent from the synthetic GWAS.

The checked-in PGEN, PVAR, and PSAM contain the same variants in the same
order. To regenerate them, run these commands in a new directory with the
fixture VCF copied there. The two imports allow PLINK to handle the
multiallelic hardcall separately from biallelic dosages.

```sh
bcftools view --no-version -m3 -Ov -o multi.vcf panel.vcf
bcftools view --no-version -M2 -Ov -o biallelic.vcf panel.vcf
plink2 --vcf multi.vcf --make-pgen --out multi
plink2 --vcf biallelic.vcf dosage=DS --hard-call-threshold 0.25 \
  --dosage-erase-threshold 0 --make-pgen --out biallelic
plink2 --pfile multi --pmerge biallelic --out panel
plink2 --pfile panel --validate
```

Generated with PLINK v2.0.0-b.1.3-devLM (26 Oct 2026) and bcftools 1.25.
Running the test suite only requires the checked-in files and bcftools.

The `multiword` fileset expands the six biallelic rows to 1,089 synthetic
samples. It exercises noncontiguous selections of 64, 65 and 129 people,
including indices on both sides of packed-word boundaries. The test compares
the PGEN results with a VCF physically subset to the same samples. To rebuild
the PGEN from its checked-in VCF, run the following command.

```sh
plink2 --vcf multiword.vcf.gz dosage=DS --hard-call-threshold 0.25 \
  --dosage-erase-threshold 0 --make-pgen --out multiword
```
