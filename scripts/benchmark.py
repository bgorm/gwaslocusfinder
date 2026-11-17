"""Measure one +gwaslocusfinder run, including its VCF/BCF and TSV writes.

Usage: python3 scripts/benchmark.py NEW_DIRECTORY +PLUGIN [general options] GWAS -- [plugin options]
The wrapper chooses the outputs and scratch directory. Inputs and flags are passed to bcftools.
"""
import csv
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time


def allocated_bytes(directory):
    total = 0
    for path in directory.rglob("*.bin"):
        try:
            total += path.stat().st_blocks * 512
        except FileNotFoundError:
            pass  # A merge pass can remove a run between listing and stat.
    return total


def main():
    if len(sys.argv) < 5 or "--" not in sys.argv[3:]:
        sys.exit(__doc__)
    destination = Path(sys.argv[1]).resolve()
    destination.mkdir()
    scratch = destination / "scratch"
    scratch.mkdir()
    command = [os.environ.get("BCFTOOLS", "bcftools"), sys.argv[2], "--no-version", "-Ob",
               "-o", destination / "index.bcf", *sys.argv[3:], "--temp-dir", scratch,
               "--leads", destination / "leads.bcf", "--independent", destination / "independent.bcf",
               "--reports-prefix", destination / "result"]
    peak_disk = 0
    with (destination / "stderr.log").open("w") as log:
        start = time.monotonic()
        with subprocess.Popen(command, stderr=log) as child:
            while True:
                peak_disk = max(peak_disk, allocated_bytes(scratch))
                try:
                    status = child.wait(timeout=0.1)
                    break
                except subprocess.TimeoutExpired:
                    pass
        wall = time.monotonic() - start
    if status:
        sys.exit(f"bcftools exited {status}; see {log.name}")
    rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    with (destination / "result.run.tsv").open() as report:
        stats = {r["PARAMETER"]: r["VALUE"] for r in csv.DictReader(report, delimiter="\t")}
    metrics = {"wall_seconds": round(wall, 3),
               "peak_rss_mib": round(rss / (1048576 if sys.platform == "darwin" else 1024), 2),
               "sampled_peak_scratch_mib": round(peak_disk / 1048576, 2)}
    for field in ("records_scanned", "reference_variants_in_windows", "reference_samples",
                  "significant_records", "loci", "genotype_reads", "tagging_records",
                  "tagging_seconds", "candidate_expansion_seconds"):
        if field in stats:
            metrics[field] = float(stats[field]) if field.endswith("_seconds") else int(stats[field])
    output = json.dumps(metrics, indent=2) + "\n"
    (destination / "metrics.json").write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
