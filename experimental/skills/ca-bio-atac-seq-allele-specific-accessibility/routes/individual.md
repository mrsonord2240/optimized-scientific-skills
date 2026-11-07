# Single-individual WASP and ASE

```bash
bash scripts/wasp_ase_pipeline.sh input.bam genotypes.vcf.gz reference.fa bowtie_index peaks.bed WASP_DIR output_dir
```

Input BAM must be paired-end, deduplicated, coordinate-sorted and indexed, with one read-group sample matching the VCF. Reference, index, VCF and peaks must share assembly and names. The pipeline requires phased heterozygotes with phase-set values for pooling, a fresh output directory, and an indexed reference; details in `references/setup-and-inputs.md`. It remaps with WASP, counts with GATK, and invokes `scripts/aggregate_peak_ase.py`. Per-site GATK itself does not require phasing, but this combined wrapper does. Do not count from the original BAM.

Default interpretation: per SNP totalCount >= 30, at least two SNPs per peak/phase set, BH p < 0.05 and absolute haplotype fraction deviation >= 0.2. Low support or smaller effects are not evidence of balance. Verify null reference fraction and report both haplotype counts and phase sets.

Done when filtered BAM and tables exist, null bias is reviewed and results are reported with limitations.
