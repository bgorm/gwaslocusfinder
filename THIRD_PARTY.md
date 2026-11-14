# Third-party sources

The project code is MIT-licensed; see LICENSE. Dependencies keep their own
licenses. The bootstrap fetches these pinned sources into vendor/:

- HTSlib 1.25: https://github.com/samtools/htslib/releases/tag/1.25
  See vendor/htslib-1.25/LICENSE and the per-directory/file notices.
- Optional pgenlib (WITH_PGEN=1) from chrchang/plink-ng commit
  4544f92e374e784a6ffd97fa88239962b81f2b62:
  https://github.com/chrchang/plink-ng/tree/4544f92e374e784a6ffd97fa88239962b81f2b62/2.0/include
  The compiled pgenlib/plink2 support files carry LGPL-3.0-or-later notices.
  See 2.0/COPYING.LESSER, 2.0/COPYING and individual source notices.

The plugin resolves HTSlib and bcftools symbols from its host and links zlib,
zstd, and system libraries. PGEN-enabled builds also include position-independent
pgenlib objects; builds with WITH_PGEN=0 neither compile nor link them.
The full pinned PLINK source archive includes SIMDe headers used on ARM; retain
their notices when redistributing those sources. GWAS expressions use the host
bcftools filter API. The plugin records and checks its bcftools/HTSlib versions.

When redistributing dependency sources or a binary containing their code, retain
the applicable licenses, notices, corresponding sources, and relinking materials.
The MIT license on project code does not replace dependency licenses.

Archive SHA-256 values:

    158c271e4e2986bb2dd11283b640edc10c56695b293df77a109ce99480303730  htslib-1.25.tar.bz2
    91bb72c289ffe820edf52eda206fd193835f2c55cbf144029207291bfb003f2c  plink-ng.tar.gz
