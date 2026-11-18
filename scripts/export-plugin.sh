#!/usr/bin/env bash
set -euo pipefail
if [[ $# != 1 || -z $1 ]]; then
    echo 'Usage: bash scripts/export-plugin.sh NEW_DIRECTORY' >&2; exit 1
fi
case "$1" in /*) destination=$1 ;; *) destination=$PWD/$1 ;; esac
if [[ -e $destination ]]; then
    echo "Destination already exists: $destination" >&2; exit 1
fi
cd "$(dirname "$0")/.."
# An explicit file list keeps notes, data, vendored checkouts, build products,
# and the internal standalone entry point out of the public source tree.
files=(
    README.md LICENSE THIRD_PARTY.md Makefile .gitignore
    src/core.hpp src/io.hpp src/storage.hpp src/reference.hpp src/entry.hpp src/engine.cpp src/plugin.cpp
    scripts/bootstrap.sh scripts/export-plugin.sh scripts/benchmark.py benchmarks/README.md
    benchmarks/measurements/provenance.json benchmarks/measurements/{84cde27,1a22388,initial}/*/metrics.json
    benchmarks/pgen-subset/{README.md,instrument.py,measure.py,pgen_probe.hpp}
    benchmarks/pgen-subset/measurements/{builds,regional-396,chr22-396,dosage-396,hardcalls-396,regional-3005,regional-3163,dosage-3005}.json
    tests/fixtures/README.md tests/fixtures/panel.vcf tests/fixtures/panel.pgen tests/fixtures/panel.pvar tests/fixtures/panel.psam
    tests/fixtures/multiword.{vcf.gz,pgen,pvar,psam}
    tests/test_plugin.py tests/test_storage.cpp .github/workflows/build.yml
)
for file in "${files[@]}"; do
    [[ -f $file ]] || { echo "Missing source file: $file" >&2; exit 1; }
done
mkdir -p "$destination"
for file in "${files[@]}"; do
    mkdir -p "$destination/$(dirname "$file")"
    cp "$file" "$destination/$file"
done
printf 'Exported plugin sources to %s\n' "$destination"
