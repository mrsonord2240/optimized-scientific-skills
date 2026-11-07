# Prepare a reference

```bash
bash scripts/prepare_reference.sh reference.fa
```

This writes a FASTA index, sequence dictionary and chromosome sizes. Input .gz must be bgzip, not plain gzip. GATK expects the dictionary named after the FASTA with .gz and the last extension removed. For assembly metadata, use `samtools dict -a GRCh38 -s "Homo sapiens" reference.fa -o reference.dict`. Compare M5 checksums and lengths for reference identity. Aligner-specific indices are handled by read-alignment/bwa-alignment (upstream only; not on this shelf) or the appropriate alignment skill.

Done when the expected files exist and a reference region can be fetched.
