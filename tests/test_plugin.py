"""Public plugin checks: requires bcftools and the plugin; the synthetic PGEN is included."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import struct
import zlib
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BCFTOOLS = os.environ.get("BCFTOOLS") or shutil.which("bcftools")
PLUGIN = Path(os.environ.get("LF_PLUGIN", ROOT / "build/gwaslocusfinder.so")).resolve()


def command(args, ok=True):
    result = subprocess.run([str(arg) for arg in args], text=True, capture_output=True)
    if (result.returncode == 0) != ok:
        raise AssertionError(f"{args}\n{result.stdout}\n{result.stderr}")
    return result


@unittest.skipUnless(BCFTOOLS and PLUGIN.is_file(), "bcftools and built plugin required")
class PublicPluginTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.supports_pgen = "--pfile <prefix>" in command([BCFTOOLS, "+" + str(PLUGIN), "--help"]).stderr

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.d = Path(self.tmp.name)
        x = [0, 0, 0, 1, 1, 1, 2, 2, 2]
        y = [0, 1, 2, 0, 1, 2, 0, 1, 2]
        z = [(a + b) / 2 for a, b in zip(x, y)]
        self.variants = [(1000, "x", x, 30, 20), (1200, "copy", x, 29, 40),
                         (2000, "z", z, 25, 5), (3000, "y", y, 20, 15)]
        self.gwas, self.ref = self.fixture("default")

    def indexed(self, path):
        output = Path(str(path) + ".bcf")
        command([BCFTOOLS, "view", "-Ob", "-o", output, "-W", path])
        return output

    def fixture(self, name, chrom="1", ref_chrom="chr1", hds=False, haploid=False):
        gwas = self.d / (name + ".gwas.vcf")
        ref = self.d / (name + ".ref.vcf")
        gwas_text = (f'##fileformat=VCFv4.2\n##contig=<ID={chrom},length=10000000>\n'
                     '##FORMAT=<ID=LP,Number=A,Type=Float,Description="Negative log10 P">\n'
                     '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n')
        ref_text = (f'##fileformat=VCFv4.2\n##contig=<ID={ref_chrom},length=10000000>\n'
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    '##FORMAT=<ID=DS,Number=A,Type=Float,Description="ALT dosage">\n'
                    '##FORMAT=<ID=HDS,Number=.,Type=Float,Description="Haplotype dosages">\n'
                    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' +
                    '\t'.join(f"s{i}" for i in range(9)) + '\n')
        for pos, ident, values, a, b in self.variants:
            gwas_text += f"{chrom}\t{pos}\t{ident}\tA\tG\t.\tPASS\t.\tLP\t{a}\t{b}\n"
            if hds:
                fmt = "HDS"
                cells = [f"{v/2},{v/2}" for v in values]
            else:
                fmt = "GT:DS"
                cells = [f"{['0/0', '0/1', '1/1'][round(v)]}:{v}" for v in values]
                if haploid:
                    cells[0] = "0:0"
            ref_text += f"{ref_chrom}\t{pos}\t{ident}\tA\tG\t.\tPASS\t.\t{fmt}\t" + '\t'.join(cells) + '\n'
        gwas.write_text(gwas_text)
        ref.write_text(ref_text)
        return self.indexed(gwas), self.indexed(ref)

    def run_plugin(self, name, general=(), specific=(), reference=True, gwas=None, ref=None, ok=True):
        prefix = self.d / name
        output = Path(str(prefix) + ".bcf")
        options = ["--reports-prefix", prefix]
        if reference:
            options += ["--vcf-ref", ref or self.ref]
        result = command([BCFTOOLS, "+" + str(PLUGIN), "-Ob", "-o", output, "-W",
                          *general, gwas or self.gwas, "--", *options, *specific], ok=ok)
        return prefix, output, result

    def ids(self, path):
        return command([BCFTOOLS, "query", "-f", "%ID\n", path]).stdout.splitlines()

    def report(self, prefix, name):
        with Path(str(prefix) + "." + name + ".tsv").open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_two_stage_clumping_and_phenotype_annotations(self):
        leads, independent = self.d / "leads.bcf", self.d / "independent.bcf"
        prefix, output, _ = self.run_plugin("loci", specific=["--leads", leads, "--independent", independent])
        self.assertEqual(self.ids(output), ["x", "copy"])
        self.assertEqual(self.ids(leads), ["x", "copy", "y"])
        self.assertEqual(self.ids(independent), ["x", "copy", "z", "y"])
        for path in [output, leads, independent]:
            self.assertTrue(Path(str(path) + ".csi").is_file())
        cells = command([BCFTOOLS, "query", "-f", "%ID[\t%LP:%LF_LOCUS]\n", output]).stdout.splitlines()
        self.assertEqual(cells, ["x\t30:1\t.:.", "copy\t.:.\t40:1"])
        self.assertEqual(len(self.report(prefix, "loci")), 2)

    def test_vcf_shards_match_single_panel(self):
        baseline, _, _ = self.run_plugin("single")
        for name, region in [("left", "chr1:1-1500"), ("right", "chr1:1501-4000")]:
            command([BCFTOOLS, "view", "-r", region, "-Ob", "-o", self.d / (name + ".bcf"), "-W", self.ref])
        shards = self.d / "refs.list"
        shards.write_text("left.bcf\nright.bcf\n")
        actual, _, _ = self.run_plugin("shards", reference=False, specific=["--vcf-ref-list", shards])
        for report in ["loci", "leads", "independent", "candidates", "excluded"]:
            self.assertEqual(self.report(baseline, report), self.report(actual, report))

    def test_filters_and_reference_sample_subset(self):
        keep = self.d / "keep.txt"
        keep.write_text(''.join(f"s{i}\n" for i in range(6)))
        prefix, output, _ = self.run_plugin("subset", general=["-s", "A", "-i", "POS>=2000"],
                                            specific=["--keep", keep, "--gwas-only"])
        self.assertEqual(self.ids(output), ["z"])
        self.assertEqual(command([BCFTOOLS, "query", "-l", output]).stdout.strip(), "A")
        run = {r["PARAMETER"]: r["VALUE"] for r in self.report(prefix, "run")}
        self.assertEqual(run["reference_samples"], "6")
        self.assertEqual(run["reference_backend"], "vcf")

    def test_hds_only_matches_ds(self):
        baseline, _, _ = self.run_plugin("ds")
        gwas, ref = self.fixture("hds", hds=True)
        actual, _, _ = self.run_plugin("hds", gwas=gwas, ref=ref)
        for report in ["loci", "leads", "independent", "candidates"]:
            self.assertEqual(self.report(baseline, report), self.report(actual, report))

    def test_chromosome_x_diploid_requirement(self):
        gwas, ref = self.fixture("x", chrom="23", ref_chrom="chrX")
        _, output, _ = self.run_plugin("x", gwas=gwas, ref=ref)
        self.assertEqual(self.ids(output), ["x", "copy"])
        gwas, ref = self.fixture("haploid", chrom="X", ref_chrom="chrX", haploid=True)
        _, _, result = self.run_plugin("haploid", gwas=gwas, ref=ref, ok=False)
        self.assertIn("requires diploid GT", result.stderr)

    def test_distance_mode_and_index_requirement(self):
        prefix, output, _ = self.run_plugin("distance", reference=False, specific=["--distance-kb", "0"])
        self.assertEqual(self.ids(output), ["x", "copy", "z", "y"])
        self.assertEqual(self.report(prefix, "leads"), [])
        _, _, result = self.run_plugin("unindexed", gwas=self.d / "default.gwas.vcf", ok=False)
        self.assertIn("index", result.stderr.lower())

    def test_vcf_forward_reads_same_position_and_cache_eviction(self):
        # More than eight metadata blocks and a 1 MiB genotype cache. The
        # same-position allele order differs from the canonical key order.
        ref = self.d / "dense.vcf"
        header = (self.d / "default.ref.vcf").read_text().split("#CHROM")[0]
        x = [0, 0, 0, 1, 1, 1, 2, 2, 2] * 12
        y = [0, 1, 2, 0, 1, 2, 0, 1, 2] * 12
        cells = lambda values: "\t".join(f"0/0:{v}" for v in values)
        with ref.open("w") as out:
            out.write(header + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
                      "\t".join(f"s{i}" for i in range(len(x))) + "\n")
            for pos in range(1000, 7000):
                if pos == 6000:
                    out.write(f"chr1\t{pos}\ty\tA\tT\t.\tPASS\t.\tGT:DS\t{cells(y)}\n")
                out.write(f"chr1\t{pos}\tx{pos}\tA\tG\t.\tPASS\t.\tGT:DS\t{cells(x)}\n")
                if pos == 6000:
                    out.write(f"chr1\t{pos}\tdiscordant\tG\tA\t.\tPASS\t.\tGT:DS\t{cells([0]*len(x))}\n")
        gwas = self.d / "dense.gwas.vcf"
        header = (self.d / "default.gwas.vcf").read_text().splitlines()
        gwas.write_text("\n".join(line for line in header if line.startswith("#")) + "\n" +
                        "1\t1000\tstart\tA\tG\t.\tPASS\t.\tLP\t18\t40\n" +
                        "1\t6000\ty\tA\tT\t.\tPASS\t.\tLP\t35\t35\n" +
                        "1\t6999\tend\tA\tG\t.\tPASS\t.\tLP\t40\t18\n")
        gwas = self.indexed(gwas)
        bcf = self.indexed(ref)
        bgzf = self.d / "dense.vcf.gz"
        command([BCFTOOLS, "view", "-Oz", "-o", bgzf, "-W", ref])
        baseline, _, _ = self.run_plugin("dense-large", gwas=gwas, ref=bcf)
        for i, panel in enumerate([bcf, bgzf]):
            actual, _, _ = self.run_plugin(f"dense-small-{i}", gwas=gwas, ref=panel,
                                          specific=["--cache-mb", "1"])
            for report in ["loci", "independent", "leads", "candidates", "excluded"]:
                self.assertEqual(self.report(actual, report), self.report(baseline, report))
            stats = {r["PARAMETER"]: r["VALUE"] for r in self.report(actual, "run")}
            self.assertEqual(stats["reference_duplicates_skipped"], "1")
            self.assertGreater(int(stats["genotype_reads"]), 6001)
            self.assertEqual(stats["reference_variants_qc_evaluated"], "6001")
            self.assertEqual(stats["reference_variants_qc_passed"], "6001")
            self.assertEqual(stats["reference_variants_qc_filtered"], "0")
            self.assertEqual(stats["reference_variants_qc_not_evaluated"], "0")
            self.assertEqual({r["N_CANDIDATE"] for r in self.report(actual, "loci")}, {"6001"})

    def test_vcf_forward_reads_do_not_cross_window_gaps(self):
        # Index a BGZF file without parsing its deliberately invalid GT fields.
        # A valid spacer allows HTSlib's usual single-position lookahead.
        header = (self.d / "default.ref.vcf").read_text().splitlines()
        text = "\n".join(line for line in header if line.startswith("#")) + "\n"
        good = "\t".join(["0/0:0", "0/1:1", "1/1:2"] * 3)
        bad = "\t".join(["invalid:0"] * 9)
        for pos, cell in [(1000, good), (1100, good), (5000, bad),
                          (20000, good), (20100, good), (1000000, bad)]:
            text += f"chr1\t{pos}\tv{pos}\tA\tG\t.\tPASS\t.\tGT:DS\t{cell}\n"
        ref = self.d / "gaps.vcf.gz"
        data = text.encode()
        with ref.open("wb") as out:
            for block in [data[i:i+32000] for i in range(0, len(data), 32000)] + [b""]:
                compressor = zlib.compressobj(wbits=-15)
                payload = compressor.compress(block) + compressor.flush()
                out.write(bytes.fromhex("1f8b08040000000000ff060042430200") +
                          struct.pack("<H", len(payload)+25) + payload +
                          struct.pack("<II", zlib.crc32(block), len(block)))
        command([BCFTOOLS, "index", "-t", ref])
        self.variants = [(1000, "first", [0]*9, 30, 30), (20000, "last", [0]*9, 20, 20)]
        gwas, _ = self.fixture("gaps-gwas")
        _, output, _ = self.run_plugin("gaps", gwas=gwas, ref=ref,
                                       specific=["--ld-kb", "0", "--merge-kb", "0"])
        self.assertEqual(self.ids(output), ["first", "last"])

    def test_reference_qc_counts_before_and_after_filters(self):
        x = [0, 1, 2] * 3
        rare = [1] + [0] * 8
        self.variants = [(1000, "common", x, 30, 30),
                         (2000, "rare", rare, 25, 25),
                         (3000, "constant", [1]*9, 24, 24),
                         (4000, "monomorphic", [0]*9, 23, 23),
                         (5000, "weak", x, 1, 1),
                         (6000, "reference_only", x, 0, 0),
                         (7000, "rare_candidate", rare, 3, 3),
                         (8000, "sparse", [0, 1, 2]*3, 22, 22)]
        self.fixture("qc")
        gwas = self.d / "qc.gwas.vcf"
        gwas.write_text("\n".join(line for line in gwas.read_text().splitlines()
                                  if line.startswith("#") or line.split("\t")[2] != "reference_only") + "\n")
        ref = self.d / "qc.ref.vcf"
        lines = ref.read_text().splitlines()
        for i, line in enumerate(lines):
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[2] == "sparse":
                fields[12:] = ["./.:."] * 6
                lines[i] = "\t".join(fields)
        ref.write_text("\n".join(lines) + "\n")
        gwas, ref = self.indexed(gwas), self.indexed(ref)
        for name, options, passed, filtered in [("unfiltered", [], 5, 2),
                                              ("maf", ["--maf", "0.1"], 3, 4),
                                              ("both", ["--maf", "0.1", "--min-call-rate", "0.9"], 2, 5)]:
            prefix, _, result = self.run_plugin("qc-"+name, gwas=gwas, ref=ref, specific=options)
            stats = {r["PARAMETER"]: r["VALUE"] for r in self.report(prefix, "run")}
            self.assertEqual(stats["reference_variants_in_windows"], "8")
            self.assertEqual(stats["reference_variants_qc_evaluated"], "7")
            self.assertEqual(stats["reference_variants_qc_passed"], str(passed))
            self.assertEqual(stats["reference_variants_qc_filtered"], str(filtered))
            self.assertEqual(stats["reference_variants_qc_not_evaluated"], "1")
            self.assertEqual(stats["genotype_reads"], "7")  # no counting-only reads
            self.assertIn("8 nearby variants before reference QC", result.stderr)
            self.assertIn(f"Reference QC: 7 variants evaluated, {passed} passed, {filtered} filtered out; 1 not evaluated (genotypes not needed)", result.stderr)

    def test_temporary_directory_cleanup_and_optional_reports(self):
        scratch = self.d / "scratch space"
        scratch.mkdir()
        sentinel = scratch / "keep.txt"
        sentinel.write_text("keep")
        _, expected, _ = self.run_plugin("scratch", specific=["--temp-dir", scratch])
        self.assertEqual(list(scratch.iterdir()), [sentinel])
        output = self.d / "without-reports.bcf"
        command([BCFTOOLS, "+" + str(PLUGIN), "-Ob", "-o", output, self.gwas, "--",
                 "--vcf-ref", self.ref, "--temp-dir", scratch])
        self.assertEqual(command([BCFTOOLS, "view", "-H", output]).stdout,
                         command([BCFTOOLS, "view", "-H", expected]).stdout)
        self.assertEqual(list(scratch.iterdir()), [sentinel])
        gwas, ref = self.fixture("bad-scratch", haploid=True)
        self.run_plugin("bad-scratch", gwas=gwas, ref=ref, specific=["--temp-dir", scratch], ok=False)
        self.assertEqual(list(scratch.iterdir()), [sentinel])
        _, _, result = self.run_plugin("missing-scratch", specific=["--temp-dir", scratch / "absent"], ok=False)
        self.assertIn("Cannot create temporary directory", result.stderr)

    def test_auxiliary_input_collisions_preserve_files(self):
        panel = ROOT / "tests/fixtures/panel"
        # Exercise primary outputs, optional streams, removed indexes, and reports.
        cases = [
            ("ids.list", "x\n", ["-i", "ID=@{file}"], [], "index"),
            ("filter samples.txt", "A\n", ["-e", "FMT/LP[@{file}:0]>100"], [], "index"),
            ("keep.txt", "s0\ns1\ns2\n", [], ["--keep", "{file}"], "index"),
            ("samples.txt", "A\n", ["-S", "{file}"], [], "index"),
            ("samples.tbi", "A\n", ["--samples-file", "^{file}"], [], "tbi"),
            ("result.run.tsv", "1\t1000\t1000\n", ["-R", "{file}"], [], "report"),
            ("regions.txt", "1\t1000\t1000\n", ["-R", "{file}"], [], "leads"),
            ("targets.txt", "1\t1000\t1000\n", ["-T", "{file}"], [], "index"),
            ("targets.csi", "1\t1500\t1600\n", ["--targets-file", "^{file}"], [], "csi"),
            ("refs.list", str(self.ref) + "\n", [], ["--vcf-ref-list", "{file}"], "leads"),
        ]
        if self.supports_pgen:
            cases.append(("pgen.list", str(panel) + "\n", [], ["--pfile-list", "{file}"], "index"))
        for name, original, general, specific, kind in cases:
            with self.subTest(file=name):
                source = self.d / name
                source.write_text(original)
                fmt = lambda args: [arg.format(file=source) for arg in args]
                prefix = self.d / "result"
                output = self.d / "index.bcf"
                output.write_bytes(b"previous index output")
                args = fmt(specific)
                if "--pfile-list" not in args and "--vcf-ref-list" not in args:
                    args += ["--vcf-ref", self.ref]
                if kind == "leads":
                    alias = self.d / (name + ".alias")
                    os.link(source, alias)
                    args += ["--leads", alias]
                elif kind in ("csi", "tbi"):
                    args += ["--independent", source.with_suffix("")]
                elif kind == "index":
                    output = source
                res = command([BCFTOOLS, "+" + str(PLUGIN), "-Ov", "-o", output, *fmt(general),
                               self.gwas, "--", *args, "--reports-prefix", prefix], ok=False)
                self.assertIn("collides with input", res.stderr)
                self.assertEqual(source.read_text(), original)
                if kind != "index":
                    self.assertEqual(output.read_bytes(), b"previous index output")

    def test_explicit_output_index_cannot_alias_input(self):
        keep = self.d / "selected.txt"
        keep.write_text("s0\ns1\ns2\n")
        before = keep.read_bytes()
        output = self.d / "main.bcf"
        output.write_bytes(b"previous variant output")
        ambiguous = f"{keep}##idx##custom.csi"
        for primary, optional in [(ambiguous, []), (output, ["--leads", ambiguous])]:
            with self.subTest(output=primary):
                result = command([BCFTOOLS, "+" + str(PLUGIN), "-Ob", "-W", "-o", primary, self.gwas,
                                  "--", "--vcf-ref", self.ref, "--keep", keep, *optional], ok=False)
                self.assertIn("##idx##", result.stderr)
                self.assertEqual(keep.read_bytes(), before)
                self.assertEqual(output.read_bytes(), b"previous variant output")

    def test_index_aliases_and_quoted_filter_literals(self):
        original_index = Path(str(self.ref) + ".csi")
        alternate = self.ref.with_suffix(".csi")
        original_index.rename(alternate)
        before = alternate.read_bytes()
        alias = self.d / "alias.bcf"
        alias.symlink_to(self.ref)
        result = command([BCFTOOLS, "+" + str(PLUGIN), "-Ov", "-o", alternate, self.gwas,
                          "--", "--vcf-ref", alias], ok=False)
        self.assertIn("collides with input", result.stderr)
        self.assertEqual(alternate.read_bytes(), before)
        # A quoted @ is part of a literal ID, not a file operand.
        output = self.d / "literal.vcf"
        output.touch()
        command([BCFTOOLS, "+" + str(PLUGIN), "-Ov", "-o", output, "-e",
                 f'ID="@{output}"', self.gwas])
        self.assertTrue(self.ids(output))

    def test_conflicting_locus_headers_fail_before_writing(self):
        template = (self.d / "default.gwas.vcf").read_text()
        prefix = self.d / "schema"
        marker = self.d / "schema.run.tsv"
        output = self.d / "schema.bcf"
        output.write_bytes(b"previous variant output")
        marker.write_bytes(b"previous completion marker")
        for tag, number, datatype in [("LF_LOCUS", "1", "String"), ("LF_START", "2", "Integer"),
                                      ("LF_END", "A", "Integer"), ("LF_INDEX_POS", ".", "Integer"),
                                      ("LF_LOCUS", "1", "Float")]:
            with self.subTest(tag=tag, number=number, datatype=datatype):
                definition = f'##FORMAT=<ID={tag},Number={number},Type={datatype},Description="Earlier definition">\n'
                raw = self.d / "conflict.vcf"
                raw.write_text(template.replace("#CHROM", definition + "#CHROM"))
                result = command([BCFTOOLS, "+" + str(PLUGIN), "-Ob", "-o", output,
                                  self.indexed(raw), "--", "--reports-prefix", prefix], ok=False)
                self.assertIn("FORMAT/" + tag, result.stderr)
                self.assertIn("Number=1,Type=Integer", result.stderr)
                self.assertEqual(output.read_bytes(), b"previous variant output")
                self.assertEqual(marker.read_bytes(), b"previous completion marker")

    def test_compatible_locus_headers_and_same_named_info(self):
        template = (self.d / "default.gwas.vcf").read_text()
        raw = self.d / "compatible.vcf"
        raw.write_text(template.replace("#CHROM", '##INFO=<ID=LF_LOCUS,Number=1,Type=String,Description="Source">\n'
                                                   '##FORMAT=<ID=LF_LOCUS,Number=1,Type=Integer,Description="Prior run">\n#CHROM')
                       .replace("\tPASS\t.\tLP\t", "\tPASS\tLF_LOCUS=original\tLP\t"))
        gwas = self.indexed(raw)
        for name in ("first", "second"):
            _, gwas, _ = self.run_plugin(name, gwas=gwas, reference=False, specific=["--distance-kb", "0"])
            query = command([BCFTOOLS, "query", "-i", "FMT/LF_LOCUS>0", "-f",
                             "%ID\t%INFO/LF_LOCUS[\t%LF_LOCUS]\n", gwas]).stdout.splitlines()
            self.assertIn("x\toriginal\t1\t1", query)
            header = command([BCFTOOLS, "view", "-h", gwas]).stdout
            self.assertEqual(header.count("##FORMAT=<ID=LF_LOCUS,"), 1)
            self.assertIn('Type=Integer,Description="Prior run"', header)

    def test_pgen_vcf_parity(self):
        if not self.supports_pgen:
            self.skipTest("PGEN decoder disabled in this build")
        panel = ROOT / "tests/fixtures/panel"
        bcf = self.indexed(shutil.copyfile(panel.with_suffix(".vcf"), self.d / "panel.vcf"))
        keep = self.d / "sample.keep"
        keep.write_text("s8\ns0\ns4\ns2\ns6\ns3\n")
        for hardcalls in (False, True):
            for subset in (False, True):
                for gwas_only in (False, True):
                    with self.subTest(hardcalls=hardcalls, subset=subset, gwas_only=gwas_only):
                        options = []
                        if hardcalls:
                            options += ["--hardcalls"]
                        if subset:
                            options += ["--keep", keep]
                        if gwas_only:
                            options += ["--gwas-only"]
                        results = []
                        for backend, path in [("--pfile", panel), ("--vcf-ref", bcf)]:
                            name = backend[2:]
                            leads, independent = self.d / (name+".leads.bcf"), self.d / (name+".ind.bcf")
                            prefix, index, _ = self.run_plugin(name, reference=False, general=["--no-version"],
                                                               specific=[backend, path, "--leads", leads,
                                                                         "--independent", independent, *options])
                            tables = {key: self.report(prefix, key) for key in
                                      ["loci", "independent", "leads", "candidates", "excluded"]}
                            streams = [command([BCFTOOLS, "view", "--no-version", "-Ov", stream]).stdout
                                       for stream in [index, leads, independent]]
                            run = {r["PARAMETER"]: r["VALUE"] for r in self.report(prefix, "run")}
                            self.assertEqual(run["reference_multiallelic_skipped"], "1")
                            self.assertEqual(run["reference_samples"], "6" if subset else "9")
                            results.append((tables, streams))
                        self.assertEqual(results[0], results[1])
                        flank = [r for r in tables["candidates"] if r["ID"] == "flank" and r["SAMPLE"] == "A"]
                        if gwas_only:
                            self.assertEqual(flank, [])
                        else:
                            hit = next(r for r in flank if r["INDEPENDENT_KEY"] == "1:1000:A:G")
                            self.assertEqual((hit["SOURCE"], hit["CALL_RATE"]), ("reference_only", "1"))
                            self.assertEqual(float(hit["R2"]) == 1, hardcalls)
                            self.assertEqual(min(int(r["START"]) for r in tables["loci"]), 700)

    def test_tagging_visits_populated_ld_windows_in_gwas_coordinates(self):
        x = [0, 1, 2] * 3
        self.variants = [(20, "distant", x, 0, 0), (900, "left", x, 2, 0),
                         (950, "filtered", x, 25, 25), (1000, "index1", x, 30, 30),
                         (1100, "right", x, 3, 0), (5000, "gap", x, 0, 0),
                         (6000, "absent", x, 30, 30), (9000, "index2", x, 28, 28)]
        self.fixture("windows", chrom="chr1", ref_chrom="1")
        ref = self.d / "windows.ref.vcf"
        ref.write_text("".join(line for line in ref.read_text().splitlines(True) if "\tabsent\t" not in line))
        gwas = self.d / "windows.gwas.vcf"
        # This deletion overlaps the first window, but its start is outside it.
        deletion = "chr1\t800\tdel\t" + "A"*151 + "\tA\t.\tPASS\t.\tLP\t4\t4\n"
        gwas.write_text(gwas.read_text().replace("chr1\t900\t", deletion + "chr1\t900\t"))
        regions = self.d / "windows.tsv"
        regions.write_text("chr1\t1000\t1000\nchr1\t6000\t6000\nchr1\t9000\t9000\n")
        prefix, _, _ = self.run_plugin("regional", gwas=self.indexed(gwas), ref=self.indexed(ref),
                                       general=["-R", regions, "-i", 'ID!="filtered"'],
                                       specific=["--ld-kb", "0.1", "--merge-kb", "0"])
        run = {r["PARAMETER"]: r["VALUE"] for r in self.report(prefix, "run")}
        self.assertEqual((run["records_scanned"], run["tagging_records"]), ("3", "5"))
        bounds = [(r["SAMPLE"], r["START"], r["END"]) for r in self.report(prefix, "loci")]
        self.assertIn(("A", "900", "1100"), bounds)
        self.assertIn(("B", "1000", "1000"), bounds)
        self.assertNotIn("filtered", {r["ID"] for r in self.report(prefix, "candidates")})

    def test_build_capabilities_and_disabled_pgen_errors(self):
        help_text = command([BCFTOOLS, "+" + str(PLUGIN), "--help"]).stderr
        enabled = "--pfile <prefix>" in help_text
        prefix, _, _ = self.run_plugin("capabilities")
        run = {r["PARAMETER"]: r["VALUE"] for r in self.report(prefix, "run")}
        self.assertEqual(run["pgen_support"], str(int(enabled)))
        if "LF_EXPECT_PGEN" in os.environ:
            self.assertEqual(enabled, os.environ["LF_EXPECT_PGEN"] == "1")
        if not enabled:
            self.assertEqual(run["pgenlib_revision"], ".")
            self.assertIn("PGEN support is disabled", help_text)
            for flag in ["--pfile", "--pfile-list"]:
                _, _, result = self.run_plugin(flag[2:], reference=False,
                                                specific=[flag, self.d / "does-not-exist"], ok=False)
                self.assertIn("PGEN support is disabled", result.stderr)
                self.assertIn("WITH_PGEN=1", result.stderr)


if __name__ == "__main__":
    unittest.main()
