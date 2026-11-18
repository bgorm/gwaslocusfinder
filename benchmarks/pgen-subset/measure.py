"""Compare release and instrumented PGEN readers on identical locus analyses."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import tempfile
import time

ROLES = ("baseline", "candidate")
OUTPUTS = ["index.bcf", "leads.bcf", "independent.bcf"] + [
    f"result.{name}.tsv" for name in ("loci", "leads", "independent", "candidates", "excluded")]


def checksum(path):
    h = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def write_sample_list(psam, destination, keep):
    with psam.open() as file:
        names = file.readline().lstrip("#").split()
        iid = names.index("IID")
        ids = [line.split()[iid] for i, line in enumerate(file) if keep(i)]
    destination.write_text("".join(s + "\n" for s in ids))
    return destination


def execute(plugin, phase, general, specific, output):
    with tempfile.TemporaryDirectory(dir=output) as tmp:
        run = Path(tmp)
        scratch = run / "scratch"
        scratch.mkdir()
        command = ["bcftools", "+" + str(plugin), "--no-version", "-Ob", "-o", run / "index.bcf",
                   *general, "--", *specific, "--temp-dir", scratch, "--leads", run / "leads.bcf",
                   "--independent", run / "independent.bcf", "--reports-prefix", run / "result"]
        env = dict(os.environ, TMPDIR=str(scratch))
        env.pop("GWASLOCUSFINDER_BENCH_VALUES", None)
        row_digest = hashlib.sha256()
        nbytes = 0
        if phase == "verify":
            env["GWASLOCUSFINDER_BENCH_VALUES"] = "1"
        with (run / "stderr.log").open("wb") as log:
            start = time.perf_counter()
            with subprocess.Popen(command, env=env, stderr=log,
                                  stdout=subprocess.PIPE if phase == "verify" else subprocess.DEVNULL) as child:
                if child.stdout is not None:
                    for block in iter(lambda: child.stdout.read(1048576), b""):
                        row_digest.update(block)
                        nbytes += len(block)
                status = child.wait()
            wall = time.perf_counter() - start
        stderr = (run / "stderr.log").read_text()
        if status:
            raise RuntimeError(stderr)
        with (run / "result.run.tsv").open() as report:
            summary = {r["PARAMETER"]: r["VALUE"] for r in csv.DictReader(report, delimiter="\t")}
        if summary["status"] != "complete" or any(scratch.iterdir()):
            raise RuntimeError("Incomplete run or scratch cleanup")
        qc = {name: int(summary[name]) for name in
              ("loci", "reference_samples", "genotype_reads", "reference_variants_qc_passed")}
        measurement = {"wall_seconds": round(wall, 6)}
        if phase != "release":
            lines = [s.split("\t") for s in stderr.splitlines() if s.startswith("PGEN_DECODE\t")]
            if len(lines) != 1 or int(lines[0][1]) != qc["genotype_reads"]:
                raise RuntimeError("Probe call count differs from the run report")
            _, _, dosage_rows, dosage_values, library, ready = lines[0]
            measurement.update(pgenlib_seconds=float(library), decoded_row_seconds=float(ready))
            if phase == "verify":
                expected_size = qc["genotype_reads"] * (4 + 4 * qc["reference_samples"])
                if nbytes != expected_size:
                    raise RuntimeError(f"Row stream contained {nbytes} bytes, expected {expected_size}")
                measurement.update(decoded_row_sha256=row_digest.hexdigest(), row_bytes=nbytes,
                                   dosage_rows=int(dosage_rows), dosage_values=int(dosage_values))
        return measurement, qc, {f: checksum(run/f) for f in OUTPUTS}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("builds", "gwas", "panel", "keep", "synthetic_panel", "synthetic_gwas", "output"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--pairs", type=int, default=6)
    args = parser.parse_args()
    if args.pairs < 2 or args.pairs % 2:
        parser.error("--pairs must be positive and even to balance the two execution orders")
    dest = args.output.resolve()
    dest.mkdir()
    plugins = {}
    commits = {}
    for role in ROLES:
        tree = (args.builds / role).resolve()
        git = ["git", "-c", f"safe.directory={tree}", "-C", tree]
        if subprocess.check_output([*git, "status", "--porcelain", "--untracked-files=no"]):
            raise RuntimeError(f"Uncommitted changes in {role}")
        commits[role] = {"commit": subprocess.check_output([*git, "rev-parse", "HEAD"], text=True).strip()}
        for phase, dirname in [("release", role), ("profile", "profile-" + role)]:
            plugin = (args.builds / dirname / "build/gwaslocusfinder.so").resolve()
            plugins[role,phase] = plugin
            commits[role][phase + "_sha256"] = checksum(plugin)
    (dest/"builds.json").write_text(json.dumps(commits, indent=2) + "\n")

    real = ["--pfile", args.panel, "--maf", "0.05", "--cache-mb", "128"]
    synth = ["--pfile", args.synthetic_panel, "--maf", "0.05", "--cache-mb", "128"]
    rkeep = ["--keep", args.keep]
    dense = write_sample_list(Path(str(args.panel)+".psam"), dest/"dense.keep", lambda i: i % 20 != 19)
    synth_psam = Path(str(args.synthetic_panel)+".psam")
    synth_keep = write_sample_list(synth_psam, dest/"synthetic.keep", lambda i: i % 8 == 0)
    synth_dense = write_sample_list(synth_psam, dest/"synthetic-dense.keep", lambda i: i % 20 != 19)
    region = ["-r", "chr22:29000001-30000000", args.gwas]
    cases = [
        ("regional-396", region, real + rkeep),
        ("chr22-396", ["-r", "chr22", args.gwas], real + rkeep),
        ("dosage-396", [args.synthetic_gwas], synth + ["--keep", synth_keep]),
        ("hardcalls-396", [args.synthetic_gwas], synth + ["--keep", synth_keep, "--hardcalls"]),
        ("regional-3005", region, real + ["--keep", dense]),
        ("regional-3163", region, real),
        ("dosage-3005", [args.synthetic_gwas], synth + ["--keep", synth_dense]),
    ]
    for name, general, specific in cases:
        print(f"Starting {name}", flush=True)
        digest = {}
        expected = None
        for role in ROLES:
            digest[role], qc, outputs = execute(plugins[role,"profile"], "verify", general, specific, dest)
            if expected is None:
                expected = qc, outputs
            elif (qc, outputs) != expected or digest[role]["decoded_row_sha256"] != digest["baseline"]["decoded_row_sha256"]:
                raise RuntimeError(f"Decoded rows, QC or outputs differ in {name}")
        result = {"qc": qc, "outputs_sha256": outputs, "row_verification": digest, "pairs": []}
        for i in range(args.pairs):
            order = list(ROLES if i % 2 == 0 else reversed(ROLES))
            pair = {"order": order}
            for phase in ("release", "profile"):
                pair[phase] = {}
                for role in order:
                    measurement, qc, outputs = execute(plugins[role,phase],phase,general,specific,dest)
                    if (qc, outputs) != expected:
                        raise RuntimeError(f"QC or selection changed in {name} pair {i+1}")
                    pair[phase][role] = measurement
            result["pairs"].append(pair)
            (dest / (name+".json")).write_text(json.dumps(result, indent=2) + "\n")
        for role in ROLES:
            wall = statistics.median(p["release"][role]["wall_seconds"] for p in result["pairs"])
            raw = statistics.median(p["profile"][role]["pgenlib_seconds"] for p in result["pairs"])
            row = statistics.median(p["profile"][role]["decoded_row_seconds"] for p in result["pairs"])
            print(f"{name} {role}: {qc['genotype_reads']} calls, wall={wall:.3f}s pgenlib={raw:.3f}s decoded={row:.3f}s",
                  flush=True)


if __name__ == "__main__":
    main()
