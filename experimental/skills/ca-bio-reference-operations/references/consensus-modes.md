# samtools consensus: Ambiguity Codes, Platform Profiles, bcftools Consensus

Moved verbatim from `SKILL.md` (2026-09-21). Read when a consensus needs IUPAC codes (`--ambig`, `--het-fract`, `--call-fract`, `--het-scale`), a platform `--config` profile, `-T` reference fill, a viral consensus command, or when deciding between `samtools consensus` (BAM) and `bcftools consensus` (reference + VCF).

### IUPAC Ambiguity for Heterozygotes
```bash
# Default Bayesian mode. --ambig is REQUIRED for IUPAC codes (R, Y, S, W, K, M, B, D, H, V);
# without it an ambiguous column is N, even when one base is 80% of the reads
samtools consensus --ambig input.bam -o consensus.fa

# Tune Bayesian het calling with --het-scale (< 1 fewer IUPAC calls, > 1 more)
samtools consensus --ambig --het-scale 0.1 input.bam -o consensus.fa

# Fixed fractional thresholds exist only in -m simple
samtools consensus -m simple --ambig --het-fract 0.2 --call-fract 0.5 input.bam -o consensus.fa
```

`--het-fract` and `--call-fract` are ignored in the default Bayesian mode (byte-identical output for 0.05, 0.9 and 0.2/0.5) and act only with `-m simple`:
- `--het-fract F`: minimum ratio of the second-most to the most common base for an IUPAC call (needs `--ambig`). Always pass it explicitly: with the flag omitted, a column with 15% minor allele stays a plain base although the help prints a default of 0.15.
- `--call-fract F`: fraction of reads that must agree on the top base, otherwise `N` (default 0.75).

`--show-ins` / `--show-del` control insertion / deletion display, not ambiguity.

### Platform-Aware Consensus
```bash
# The default Bayesian algorithm needs no --config; platform-specific profiles (samtools 1.17+; list them with `samtools help consensus`)
samtools consensus --config hifi       input.bam -o consensus.fa   # PacBio HiFi
samtools consensus --config r10.4_sup  input.bam -o consensus.fa   # ONT R10.4+ (r10.4_dup for duplex)
samtools consensus --config ultima     input.bam -o consensus.fa   # Ultima Genomics
samtools consensus --config hiseq      input.bam -o consensus.fa   # Illumina

# Report the ref base at columns with no reads (depth 0; -T added in samtools 1.22; bases keep the FASTA's case)
samtools consensus -T ref.fa input.bam -o consensus.fa
```
`-T` fills only depth-0 columns: columns below `-d` and ambiguous columns stay `N` (planted 450-`N` result at `-d 3`: 300 zero-depth columns filled, 150 below `-d` stay `N`).

### samtools consensus vs bcftools consensus

Different operations -- conflating them produces nonsense:

| Tool | Input | Output | Use case |
|------|-------|--------|----------|
| `samtools consensus` | BAM | Consensus FASTA derived from reads (Bayesian) | Viral, de novo / amplicon, low-coverage species |
| `bcftools consensus` | reference + VCF | Reference with VCF variants applied | Apply called variants (haplotype reconstruction, custom ref for re-mapping) |

For viral consensus from BAM:
```bash
# Modern: samtools consensus
# (--show-del yes would write '*' into the FASTA, so the default no is kept)
samtools consensus --config hiseq -d 10 --ambig -a input.bam -o consensus.fa

# Apply called variants to reference (different question)
bcftools consensus -f reference.fa variants.vcf.gz -o sample_consensus.fa
bcftools consensus -f reference.fa -H 1 phased.vcf.gz -o haplotype1.fa   # phased haplotype 1
```

With a genotyped (FORMAT/GT) VCF and no `-H`, bcftools 1.24 writes heterozygous SNPs as IUPAC codes; use `-H 1` / `-H 2` for one haplotype, or `-H A` (or `-s -`) to apply every ALT allele.

`samtools consensus` is not iterative and is not an assembly-polishing tool.
