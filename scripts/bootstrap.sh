#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
with_pgen=${WITH_PGEN:-1}
for arg in "$@"; do
    case "$arg" in
        --without-pgen) with_pgen=0 ;;
        --plugin) ;; # accepted for compatibility; the plugin is now the default
        -h|--help) echo 'Usage: bash scripts/bootstrap.sh [--without-pgen]'; exit 0 ;;
        *) echo "Unknown option: $arg" >&2; exit 1 ;;
    esac
done
case "$with_pgen" in 0|1) ;; *) echo 'WITH_PGEN must be 0 or 1' >&2; exit 1 ;; esac
verify_sha256() {
    local expected=$1 archive=$2 actual
    if command -v sha256sum >/dev/null 2>&1; then
        actual=$(sha256sum "$archive")
    elif command -v shasum >/dev/null 2>&1; then
        actual=$(shasum -a 256 "$archive")
    else
        echo 'A SHA-256 checker is required (sha256sum or shasum)' >&2; exit 1
    fi
    if [[ ${actual%% *} != "$expected" ]]; then
        echo "SHA-256 mismatch: $archive" >&2; exit 1
    fi
}
htslib_dir=${PLUGIN_HTSLIB:-${HTSLIB:-vendor/htslib-1.25}}
if [[ ! -f $htslib_dir/htslib/hts.h ]]; then
    if [[ $htslib_dir != vendor/htslib-1.25 ]]; then
        echo "HTSlib headers not found in $htslib_dir" >&2; exit 1
    fi
    mkdir -p vendor
    curl -fL --retry 2 https://github.com/samtools/htslib/releases/download/1.25/htslib-1.25.tar.bz2 -o vendor/htslib-1.25.tar.bz2
    verify_sha256 158c271e4e2986bb2dd11283b640edc10c56695b293df77a109ce99480303730 vendor/htslib-1.25.tar.bz2
    tar -xjf vendor/htslib-1.25.tar.bz2 -C vendor
fi
revision=4544f92e374e784a6ffd97fa88239962b81f2b62
pgen_dir=${PGENLIB:-vendor/plink-ng-$revision/2.0/include}
if [[ $with_pgen == 1 && ! -f $pgen_dir/pgenlib_read.h ]]; then
    if [[ $pgen_dir != vendor/plink-ng-$revision/2.0/include ]]; then
        echo "pgenlib sources not found in $pgen_dir" >&2; exit 1
    fi
    mkdir -p vendor
    curl -fL --retry 2 "https://codeload.github.com/chrchang/plink-ng/tar.gz/$revision" -o vendor/plink-ng.tar.gz
    verify_sha256 91bb72c289ffe820edf52eda206fd193835f2c55cbf144029207291bfb003f2c vendor/plink-ng.tar.gz
    tar -xzf vendor/plink-ng.tar.gz -C vendor
fi
# The plugin needs HTSlib headers; it uses the host's HTSlib at runtime.
make -j"${JOBS:-4}" plugin WITH_PGEN="$with_pgen" PLUGIN_HTSLIB="$htslib_dir" PGENLIB="$pgen_dir"
