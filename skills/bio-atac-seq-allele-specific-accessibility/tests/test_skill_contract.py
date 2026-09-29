import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class SkillContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pipeline = (ROOT / "scripts" / "wasp_ase_pipeline.sh").read_text(encoding="utf-8")
        cls.skill = (ROOT / "SKILL.md").read_text(encoding="utf-8")
        cls.methods = (ROOT / "references" / "method-selection-and-failures.md").read_text(encoding="utf-8")
        cls.provenance = (ROOT / "references" / "provenance.md").read_text(encoding="utf-8")

    def test_final_bam_is_sorted_read_grouped_indexed_and_nonempty(self):
        final_merge = self.pipeline.index('samtools merge -f "$STAGE/wasp/${SAMPLE}.merged.unsorted.bam"')
        operations = [
            "samtools sort -o", "samtools addreplacerg", "samtools index",
            "GATK ASEReadCounter wrote no data rows",
        ]
        positions = [self.pipeline.index(operation, final_merge) for operation in operations]
        self.assertEqual(positions, sorted(positions))
        self.assertIn('FINAL_SORT_ORDER', self.pipeline)
        self.assertIn('FINAL_SAMPLES', self.pipeline)

    def test_bam_sample_is_resolved_before_target_vcf_filter(self):
        resolve = self.pipeline.index("BAM_SAMPLES")
        select = self.pipeline.index('bcftools view --samples "$SAMPLE"')
        target_output = self.pipeline.index('-o "$TARGET_VCF"', select)
        heterozygous = self.pipeline.index("-i 'GT=\"het\"'", target_output)
        self.assertLess(resolve, select)
        self.assertLess(target_output, heterozygous)
        self.assertIn('-o "$FILTERED_VCF" "$TARGET_VCF"', self.pipeline)
        self.assertIn("must match exactly one VCF sample", self.pipeline)

    def test_run_is_staged_and_existing_output_is_refused(self):
        self.assertIn('if [[ -e "$OUTDIR" ]]', self.pipeline)
        self.assertIn('STAGE=$(mktemp -d "${OUTDIR}.tmp.XXXXXX")', self.pipeline)
        self.assertIn('mv -- "$STAGE" "$OUTDIR"', self.pipeline)
        unquoted_redirects = re.findall(r">\s+\$(?:OUTDIR|BAM|VCF|GENOME_FA|PEAKS|WASP_DIR)\b", self.pipeline)
        self.assertEqual(unquoted_redirects, [])

    def test_matrixeqtl_and_quasar_are_explicit_routing_only(self):
        combined = self.skill + self.methods
        self.assertIn("method-selection routes only", combined)
        self.assertIn("no executable QuASAR", combined)
        self.assertIn("no executable MatrixEQTL", combined)

    def test_scientific_claims_are_conditioned(self):
        combined = self.skill + self.methods
        self.assertNotIn("1.5-3x", combined)
        self.assertNotIn("above 70%", combined)
        self.assertNotIn("p < 1e-5", combined)
        self.assertIn("study-specific FDR or permutation", combined)

    def test_provenance_uses_canonical_git_blob_ids(self):
        for blob in (
            "e0941c7bde0a2e4540e6f1b962575bfcfa77254a",
            "f3e7ccbebed3654a2efa24b5a62b58e43719c0cc",
            "e8920f67859dede6e3f9b49793047a7a153e91da",
        ):
            self.assertIn(blob, self.provenance)
        self.assertIn("may reflect its newline normalization", self.provenance)


if __name__ == "__main__":
    unittest.main()
