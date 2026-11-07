# Consensus and comparison

```bash
samtools consensus input.bam -o consensus.fa
```

Use `-r chr1:1000-2000` for a region, `-d 5` for a minimum depth of five, and `-f fastq -l 0` for unwrapped FASTQ. `-a` pads covered contigs to reference ends, `-aa` also emits empty contigs. For reference-aligned columns use `-a --show-del yes --show-ins no`. Choose IUPAC, platform profiles, reference fill, viral consensus, or a VCF-based consensus with `references/consensus-modes.md`.

To compare a consensus to reference, `minimap2 -a reference.fa consensus.fa > comparison.sam`. The teaching-only Python majority-vote method is not quality-aware: `python scripts/pysam_consensus.py compare input.bam reference.fa chr22 1951 4617 --min-depth 3`; its interval is 0-based half-open. See `references/python-consensus.md`.

Done when the requested consensus or comparison exists and its coverage and ambiguity choices are reported.
