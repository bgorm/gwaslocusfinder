"""Public plugin checks: no standalone executable, PLINK, or PGEN files required."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
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

    def test_build_capabilities_and_disabled_pgen_errors(self):
        help_text = command([BCFTOOLS, "+" + str(PLUGIN), "--help"]).stderr
        enabled = "--pfile <prefix>" in help_text
        prefix, _, _ = self.run_plugin("capabilities")
        run = {r["PARAMETER"]: r["VALUE"] for r in self.report(prefix, "run")}
        self.assertEqual(run["pgen_support"], str(int(enabled)))
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
