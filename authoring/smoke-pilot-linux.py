"""Offline BAM/VCF fixture checks in WSL; accepts an exported snapshot directory.

Export with smoke-pilot.mjs --export-only and the four IDs below. This runner
checks exported SHA-256 values before executing any bundled code. It changes
only the uniquely generated scratch directory, and preserves outputs for review.
"""
import hashlib
import json
import pathlib
import subprocess
import sys

import pysam

root = pathlib.Path(sys.argv[1]).resolve()
manifest = json.loads((root / "source-manifest.json").read_text())
for name, expected in manifest["files"].items():
    actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
    if actual != expected:
        raise ValueError(f"Export hash mismatch: {name}")
tests = []


def run(skill, script, args, check):
    directory = root / skill
    result = subprocess.run(args, cwd=directory, capture_output=True, text=True, timeout=60)
    detail = None
    try:
        assert result.returncode == 0, f"exit {result.returncode}"
        check(directory, result.stdout)
        status = "passed"
    except (AssertionError, OSError, ValueError) as error:
        status, detail = "failed", str(error)
    tests.append({"id": skill, "script": script, "args": args, "status": status,
                  "exit_code": result.returncode, "error": detail,
                  "source_sha256": manifest["files"][f"{skill}/{script}"],
                  "stdout": result.stdout, "stderr": result.stderr})


def count_bam(file):
    with pysam.AlignmentFile(file, "rb") as bam:
        return sum(1 for _ in bam.fetch(until_eof=True))


def bam_fixture(directory):
    file = directory / "smoke.bam"
    with pysam.AlignmentFile(file, "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"},
                                               "SQ": [{"SN": "chr1", "LN": 2000}]}) as bam:
        for i in range(200):
            read = pysam.AlignedSegment()
            read.query_name, read.query_sequence = f"read{i}", "A" * 10
            read.flag, read.reference_id, read.reference_start = 16 if i % 2 else 0, 0, i * 5
            read.mapping_quality, read.cigarstring = 10 if i % 5 == 0 else 60, "10M"
            read.query_qualities = pysam.qualitystring_to_array("I" * 10)
            bam.write(read)
    pysam.index(str(file))


def expect_bam(file, count):
    def check(directory, stdout):
        assert count_bam(directory / file) == count, f"Expected {count} reads"
    return check


filter_id = "bio-alignment-filtering"
bam_fixture(root / filter_id)
(root / filter_id / "smoke.bed").write_text("chr1\t0\t60\nchr1\t50\t100\n")
run(filter_id, "examples/filter_bam.py",
    [sys.executable, "-B", "examples/filter_bam.py", "smoke.bam", "filtered.bam", "-q", "30"],
    expect_bam("filtered.bam", 160))
run(filter_id, "scripts/filter_by_bed.py",
    [sys.executable, "-B", "scripts/filter_by_bed.py", "smoke.bam", "smoke.bed", "regions.bam"],
    expect_bam("regions.bam", 20))

validation_id = "bio-alignment-validation"
bam_fixture(root / validation_id)


def check_validation(directory, stdout):
    assert "Mapped: 200 (100.00%)" in stdout
    assert "Forward fraction F/(F+R): 0.500" in stdout
    assert "All metrics within normal range" in stdout


run(validation_id, "examples/validate_alignment.py",
    [sys.executable, "-B", "examples/validate_alignment.py", "smoke.bam"], check_validation)

normalization_id = "bio-variant-normalization"
directory = root / normalization_id
(directory / "reference.fa").write_text(">chr1\n" + "A" * 1200 + "\n")
pysam.faidx(str(directory / "reference.fa"))
(directory / "input.vcf").write_text(
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1200>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr1\t100\t.\tA\tC,G\t60\tPASS\t.\nchr1\t200\t.\tAA\tCG\t60\tPASS\t.\n")


def check_normalized(directory, stdout):
    with pysam.VariantFile(directory / "normalized.vcf.gz") as variants:
        records = list(variants.fetch())
    assert len(records) == 4
    assert all(len(record.ref) == 1 and len(record.alts) == 1 and len(record.alts[0]) == 1 for record in records)


run(normalization_id, "examples/normalize_vcf.sh",
    ["bash", "examples/normalize_vcf.sh", "reference.fa", "input.vcf", "normalized.vcf.gz"], check_normalized)

variant_filter_id = "bio-variant-calling-filtering-best-practices"
directory = root / variant_filter_id
info_fields = ["QD", "FS", "MQ", "MQRankSum", "ReadPosRankSum", "SOR"]
(directory / "input.vcf").write_text(
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1200>\n" +
    "".join(f'##INFO=<ID={name},Number=1,Type=Float,Description="Fixture {name}">\n' for name in info_fields) +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t100\tkeep_missing\tA\tG\t60\tPASS\tQD=10;FS=1;MQ=60;SOR=1\n" +
    "chr1\t110\tdrop_lowqual\tA\tC\t10\tPASS\tQD=10;FS=1;MQ=60;SOR=1\n" +
    "chr1\t120\tkeep_indel\tAA\tA\t60\tPASS\tQD=10;FS=100;SOR=1\n")


def check_filter(directory, stdout):
    with pysam.VariantFile(directory / "filtered_all_filtered.vcf.gz") as variants:
        ids = {record.id for record in variants.fetch()}
    assert ids == {"keep_missing", "keep_indel"}, str(ids)


run(variant_filter_id, "examples/filter_variants.sh",
    ["bash", "examples/filter_variants.sh", "input.vcf", "filtered"], check_filter)
result = {"schema_version": 1, "scope": "provider_offline_smoke_only_not_aipoch_approval",
          "source_commit": manifest["source_commit"], "python_version": sys.version,
          "pysam_version": pysam.__version__, "scratch_directory": str(root),
          "samtools_version": subprocess.check_output(["samtools", "--version"], text=True, errors="replace").splitlines()[0],
          "bcftools_version": subprocess.check_output(["bcftools", "--version"], text=True).splitlines()[0],
          "tests": tests, "passed": sum(test["status"] == "passed" for test in tests),
          "failed": sum(test["status"] != "passed" for test in tests)}
print(json.dumps(result, indent=2))
sys.exit(bool(result["failed"]))
