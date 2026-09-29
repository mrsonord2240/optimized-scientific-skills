#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


def hyb_row(read_id: str, mirna_first: bool = False, target: str = "ENST1_gene_mRNA") -> str:
    mirna = ["MIMAT1_MirBase_miR-1_microRNA", "1", "22", "1", "22", "1e-4"]
    mrna = [target, "23", "40", "101", "118", "2e-4"]
    fields = [read_id, "ACGT" * 10, "."] + (mirna + mrna if mirna_first else mrna + mirna) + [""]
    return "\t".join(fields) + "\n"


class ScriptTests(unittest.TestCase):
    def run_cmd(self, *args: object, check: bool = True, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
        return subprocess.run([str(arg) for arg in args], text=True, capture_output=True, check=check, env=env)

    def test_consensus_orientation_expression_and_ambiguity(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); run1 = root / "run1.hyb"; run2 = root / "run2.hyb"
            run1.write_text(hyb_row("stable", False) + hyb_row("unstable", True, "ENST2_gene_mRNA"), encoding="utf-8")
            run2.write_text(hyb_row("stable", False) + hyb_row("unstable", True, "ENST3_gene_mRNA"), encoding="utf-8")
            expression = root / "expression.tsv"
            expression.write_text("mirna_id\texpression_value\texpression_unit\texpression_source\nMIMAT1_MirBase_miR-1_microRNA\t150\tTPM\tmatched-small-rna\n", encoding="utf-8")
            sites, targets, excluded, support, manifest = (root / name for name in ("sites.tsv", "targets.tsv", "excluded.tsv", "support.tsv", "manifest.json"))
            self.run_cmd("python3", SCRIPTS / "consensus_hyb.py", "--hyb", run1, run2,
                         "--sites", sites, "--targets", targets, "--excluded", excluded,
                         "--support", support, "--manifest", manifest, "--expression", expression,
                         "--expression-threshold", "100", "--hyb-commit", "abc",
                         "--hyb-db", "db", "--run-id", "test", "--reads", run1)
            with sites.open(encoding="utf-8") as handle:
                site_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["read_id"] for row in site_rows], ["stable"])
            self.assertEqual(site_rows[0]["source_orientation"], "target-first")
            self.assertEqual(site_rows[0]["expression_source"], "matched-small-rna")
            self.assertIn("unstable_assignment", excluded.read_text(encoding="utf-8"))
            with support.open(encoding="utf-8") as handle:
                support_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(support_rows), 3)
            self.assertEqual([row["runs_supporting_assignment"] for row in support_rows], ["2", "1", "1"])
            self.assertEqual([row["accepted"] for row in support_rows], ["true", "false", "false"])
            self.assertEqual(json.loads(manifest.read_text())["accepted_rows"], 1)

    def test_consensus_rejects_wrong_schema(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); bad = root / "bad.hyb"; bad.write_text("a\tb\tc\n", encoding="utf-8")
            result = self.run_cmd("python3", SCRIPTS / "consensus_hyb.py", "--hyb", bad, bad,
                                  "--sites", root / "s", "--targets", root / "t", "--excluded", root / "e",
                                  "--support", root / "u", "--manifest", root / "m", "--hyb-commit", "abc", "--hyb-db", "db",
                                  "--run-id", "test", "--reads", bad, check=False)
            self.assertNotEqual(result.returncode, 0); self.assertIn("expected 16", result.stderr)

    def test_expression_values_must_be_finite_and_numeric(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); run1 = root / "run1.hyb"; run2 = root / "run2.hyb"
            high_mirna = "MIMAT1_MirBase_miR-1_microRNA"
            low_mirna = "MIMAT2_MirBase_miR-2_microRNA"

            def row(read_id: str, mirna: str) -> str:
                mrna = ["ENST1_gene_mRNA", "23", "40", "101", "118", "2e-4"]
                mirna_segment = [mirna, "1", "22", "1", "22", "1e-4"]
                return "\t".join([read_id, "ACGT" * 10, ".", *mrna, *mirna_segment, ""]) + "\n"

            fixture_rows = row("above-threshold", high_mirna) + row("below-threshold", low_mirna)
            run1.write_text(fixture_rows, encoding="utf-8"); run2.write_text(fixture_rows, encoding="utf-8")
            expression = root / "expression.tsv"
            expression.write_text(
                "mirna_id\texpression_value\texpression_unit\texpression_source\n"
                f"{high_mirna}\t150\tTPM\tmatched-small-rna\n"
                f"{low_mirna}\t50\tTPM\tmatched-small-rna\n", encoding="utf-8")
            outputs = {key: root / value for key, value in {
                "sites": "sites.tsv", "targets": "targets.tsv", "excluded": "excluded.tsv",
                "support": "support.tsv", "manifest": "manifest.json",
            }.items()}
            self.run_cmd("python3", SCRIPTS / "consensus_hyb.py", "--hyb", run1, run2,
                         "--sites", outputs["sites"], "--targets", outputs["targets"],
                         "--excluded", outputs["excluded"], "--support", outputs["support"],
                         "--manifest", outputs["manifest"], "--reads", run1,
                         "--expression", expression, "--expression-threshold", "100",
                         "--hyb-commit", "abc", "--hyb-db", "db", "--run-id", "finite-expression")
            with outputs["sites"].open(encoding="utf-8") as handle:
                sites = list(csv.DictReader(handle, delimiter="\t"))
            with outputs["excluded"].open(encoding="utf-8") as handle:
                excluded = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([site["read_id"] for site in sites], ["above-threshold"])
            self.assertEqual(sites[0]["expression_value"], "150.0")
            self.assertEqual(excluded, [{"read_id": "below-threshold", "reason": "below_expression_threshold",
                                         "runs_present": "2", "runs_required": "2"}])
            self.assertEqual(json.loads(outputs["manifest"].read_text(encoding="utf-8"))["accepted_rows"], 1)

            invalid_values = ("NaN", "nan", "+Inf", "Infinity", "nAn", "-inf", "-Infinity")
            for index, value in enumerate(invalid_values):
                with self.subTest(expression_value=value):
                    invalid_root = root / f"nonfinite-{index}"; invalid_root.mkdir()
                    invalid_expression = invalid_root / "expression.tsv"
                    invalid_expression.write_text(
                        "mirna_id\texpression_value\texpression_unit\texpression_source\n"
                        f"{high_mirna}\t{value}\tTPM\tmatched-small-rna\n", encoding="utf-8")
                    invalid_outputs = {key: invalid_root / path.name for key, path in outputs.items()}
                    result = self.run_cmd(
                        "python3", SCRIPTS / "consensus_hyb.py", "--hyb", run1, run2,
                        "--sites", invalid_outputs["sites"], "--targets", invalid_outputs["targets"],
                        "--excluded", invalid_outputs["excluded"], "--support", invalid_outputs["support"],
                        "--manifest", invalid_outputs["manifest"], "--reads", run1,
                        "--expression", invalid_expression, "--expression-threshold", "100",
                        "--hyb-commit", "abc", "--hyb-db", "db", "--run-id", "invalid-expression",
                        check=False)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(f"{invalid_expression}:2: expression_value must be finite", result.stderr)
                    self.assertTrue(all(not path.exists() for path in invalid_outputs.values()))

            for index, (value, header) in enumerate((("not-a-number", "normal"), ("", "normal"), ("", "short"))):
                with self.subTest(malformed_expression_value=value, header=header):
                    invalid_root = root / f"malformed-{index}"; invalid_root.mkdir()
                    invalid_expression = invalid_root / "expression.tsv"
                    expression_header = (
                        "mirna_id\texpression_value\texpression_unit\texpression_source\n"
                        if header == "normal" else
                        "mirna_id\texpression_value\texpression_unit\texpression_source\n"
                    )
                    value_row = (
                        f"{high_mirna}\t{value}\tTPM\tmatched-small-rna\n"
                        if header == "normal" else f"{high_mirna}\n"
                    )
                    invalid_expression.write_text(expression_header + value_row, encoding="utf-8")
                    invalid_outputs = {key: invalid_root / path.name for key, path in outputs.items()}
                    result = self.run_cmd(
                        "python3", SCRIPTS / "consensus_hyb.py", "--hyb", run1, run2,
                        "--sites", invalid_outputs["sites"], "--targets", invalid_outputs["targets"],
                        "--excluded", invalid_outputs["excluded"], "--support", invalid_outputs["support"],
                        "--manifest", invalid_outputs["manifest"], "--reads", run1,
                        "--expression", invalid_expression, "--expression-threshold", "100",
                        "--hyb-commit", "abc", "--hyb-db", "db", "--run-id", "invalid-expression",
                        check=False)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(f"{invalid_expression}:2: invalid expression_value", result.stderr)
                    self.assertTrue(all(not path.exists() for path in invalid_outputs.values()))

            bad_threshold_root = root / "bad-threshold"; bad_threshold_root.mkdir()
            bad_threshold_outputs = {key: bad_threshold_root / path.name for key, path in outputs.items()}
            result = self.run_cmd(
                "python3", SCRIPTS / "consensus_hyb.py", "--hyb", run1, run2,
                "--sites", bad_threshold_outputs["sites"], "--targets", bad_threshold_outputs["targets"],
                "--excluded", bad_threshold_outputs["excluded"], "--support", bad_threshold_outputs["support"],
                "--manifest", bad_threshold_outputs["manifest"], "--reads", run1,
                "--expression", expression, "--expression-threshold", "NaN",
                "--hyb-commit", "abc", "--hyb-db", "db", "--run-id", "invalid-threshold", check=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("expression threshold must be finite", result.stderr)
            self.assertTrue(all(not path.exists() for path in bad_threshold_outputs.values()))

    def test_targetscan_plus_minus_and_release_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); sites = root / "sites.tsv"; maps = root / "maps.tsv"
            sites.write_text(
                "transcript_id\tutr_start_1\tutr_end_1\tmirna_family\tsite_type\n"
                "txPlus\t8\t14\tmiR-1\t8mer\n"
                "txMinus\t8\t14\tmiR-1\t8mer\n", encoding="utf-8")
            maps.write_text(
                "transcript_id\tchrom\tstrand\tutr_exon_starts_0\tutr_exon_ends_0\tassembly\tannotation_release\ttargetscan_release\n"
                "txPlus\tchr1\t+\t100,200\t110,210\tGRCh38\tGENCODEv49\t8.0\n"
                "txMinus\tchr2\t-\t300,400\t310,410\tGRCh38\tGENCODEv49\t8.0\n", encoding="utf-8")
            bed, manifest = root / "sites.bed", root / "manifest.json"
            self.run_cmd("python3", SCRIPTS / "targetscan_sites_to_bed12.py", "--sites", sites,
                         "--utr-map", maps, "--output", bed, "--manifest", manifest,
                         "--assembly", "GRCh38", "--annotation-release", "GENCODEv49", "--targetscan-release", "8.0")
            rows = [line.split("\t") for line in bed.read_text().splitlines()]
            self.assertEqual(rows[0][5], "+"); self.assertEqual(rows[0][9], "2")
            self.assertEqual(rows[1][5], "-"); self.assertEqual(rows[1][9], "2")
            self.assertEqual(json.loads(manifest.read_text())["strand_preserved"], True)

    def test_wrapper_quotes_paths_and_publishes_atomically(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); hyb_home = root / "hyb home"; db = hyb_home / "data" / "db"
            db.mkdir(parents=True); (db / "testdb.1.bt2").write_text("index")
            reads = root / "reads * literal.fastq"; reads.write_text("@r\nACGT\n+\nIIII\n")
            fake = root / "fake hyb.sh"
            fake.write_text("#!/usr/bin/env bash\nset -euo pipefail\n[[ ${PERL_HASH_SEED-} == 0 && ${PERL_PERTURB_KEYS-} == 0 && ${PYTHONHASHSEED-} == 0 && ${LC_ALL-} == C && ${TZ-} == UTC ]] || exit 44\nid=; db=\nfor arg in \"$@\"; do case \"$arg\" in id=*) id=${arg#id=};; db=*) db=${arg#db=};; esac; done\nprintf 'r1\\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\t.\\tENST1_gene_mRNA\\t23\\t40\\t101\\t118\\t2e-4\\tMIMAT1_MirBase_miR-1_microRNA\\t1\\t22\\t1\\t22\\t1e-4\\t\\n' > \"${id}_comp_${db}_hybrids_ua.hyb\"\n", encoding="utf-8")
            fake.chmod(0o755); output = root / "output with spaces"
            env = os.environ.copy(); env["HYB_HOME"] = str(hyb_home)
            self.run_cmd("bash", SCRIPTS / "run_chimeric_eclip.sh", "--reads", reads, "--hyb-db", "testdb",
                         "--run-id", "sample", "--output-dir", output, "--hyb-bin", fake, env=env)
            self.assertTrue((output / "manifest.json").is_file())
            self.assertEqual(json.loads((output / "manifest.json").read_text())["accepted_rows"], 1)
            self.assertTrue((output / "support.tsv").is_file())
            repeated = root / "repeated output"
            self.run_cmd("bash", SCRIPTS / "run_chimeric_eclip.sh", "--reads", reads, "--hyb-db", "testdb",
                         "--run-id", "sample", "--output-dir", repeated, "--hyb-bin", fake, env=env)
            for name in ("sites.tsv", "targets.tsv", "excluded.tsv", "support.tsv", "manifest.json"):
                self.assertEqual((output / name).read_bytes(), (repeated / name).read_bytes())
            second = self.run_cmd("bash", SCRIPTS / "run_chimeric_eclip.sh", "--reads", reads, "--hyb-db", "testdb",
                                  "--run-id", "sample", "--output-dir", output, "--hyb-bin", fake,
                                  check=False, env=env)
            self.assertEqual(second.returncode, 73); self.assertTrue((output / "manifest.json").is_file())

    def test_wrapper_propagates_hyb_status_and_preserves_existing_output(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); hyb_home = root / "hyb"; db = hyb_home / "data" / "db"
            db.mkdir(parents=True); (db / "testdb.1.bt2").write_text("index")
            reads = root / "reads.fastq"; reads.write_text("@r\nACGT\n+\nIIII\n")
            fake = root / "failing-hyb"
            fake.write_text("#!/usr/bin/env bash\necho 'controlled Hyb failure' >&2\nexit 42\n", encoding="utf-8"); fake.chmod(0o755)
            output = root / "output"; output.mkdir(); (output / "sentinel.txt").write_text("preserve")
            env = os.environ.copy(); env["HYB_HOME"] = str(hyb_home)
            result = self.run_cmd("bash", SCRIPTS / "run_chimeric_eclip.sh", "--reads", reads,
                                  "--hyb-db", "testdb", "--run-id", "sample", "--output-dir", output,
                                  "--hyb-bin", fake, "--replace", check=False, env=env)
            self.assertEqual(result.returncode, 42); self.assertIn("controlled Hyb failure", result.stderr)
            self.assertEqual((output / "sentinel.txt").read_text(), "preserve")

    def test_targeted_umi_requires_declared_length_and_records_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp); r1 = root / "R1.fastq.gz"; r2 = root / "R2.fastq.gz"
            with gzip.open(r1, "wt", encoding="utf-8") as handle:
                handle.write("@read1 description\nAACCGGTT\n+\nIIIIIIII\n")
            with gzip.open(r2, "wt", encoding="utf-8") as handle:
                handle.write("@read1 description\nACGTACGTACGG\n+\nIIIIIIIIIIII\n")

            missing = self.run_cmd("python3", SCRIPTS / "extract_targeted_umi.py",
                                   "--read1", r1, "--read2", r2,
                                   "--output-fastq", root / "missing.fastq", "--manifest", root / "missing.json",
                                   "--library-id", "fixture", "--protocol-source", "fixture-protocol-v1",
                                   check=False)
            self.assertNotEqual(missing.returncode, 0); self.assertIn("--umi-length", missing.stderr)

            invalid_cases = (
                ("zero", "0", "fixture", "fixture-protocol-v1", "between 1 and 64"),
                ("too-long", "65", "fixture", "fixture-protocol-v1", "between 1 and 64"),
                ("not-integer", "nine", "fixture", "fixture-protocol-v1", "invalid int value"),
                ("blank-library", "9", "   ", "fixture-protocol-v1", "must be non-empty"),
                ("blank-protocol", "9", "fixture", "   ", "must be non-empty"),
            )
            for label, length, library_id, protocol_source, diagnostic in invalid_cases:
                with self.subTest(label=label):
                    output = root / f"invalid-{label}.fastq"
                    manifest = root / f"invalid-{label}.json"
                    result = self.run_cmd(
                        "python3", SCRIPTS / "extract_targeted_umi.py",
                        "--read1", r1, "--read2", r2,
                        "--output-fastq", output, "--manifest", manifest,
                        "--umi-length", length, "--library-id", library_id,
                        "--protocol-source", protocol_source, check=False,
                    )
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(diagnostic, result.stderr)
                    self.assertFalse(output.exists())
                    self.assertFalse(manifest.exists())

            for length, expected in ((9, "ACGTACGTA"), (10, "ACGTACGTAC")):
                output = root / f"umi-{length}.fastq"; manifest = root / f"umi-{length}.json"
                self.run_cmd("python3", SCRIPTS / "extract_targeted_umi.py",
                             "--read1", r1, "--read2", r2, "--output-fastq", output,
                             "--manifest", manifest, "--umi-length", length,
                             "--library-id", f"fixture-{length}",
                             "--protocol-source", f"fixture-protocol-{length}nt")
                self.assertTrue(output.read_text(encoding="utf-8").splitlines()[0].startswith(f"@read1_{expected} "))
                record = json.loads(manifest.read_text(encoding="utf-8"))
                self.assertEqual(record["declared_umi_length"], length)
                self.assertEqual(record["protocol_source"], f"fixture-protocol-{length}nt")
                self.assertEqual(record["output"]["reads"], 1)


if __name__ == "__main__": unittest.main(verbosity=2)
