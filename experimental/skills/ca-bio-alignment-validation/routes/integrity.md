# File and reference integrity

```bash
samtools quickcheck -v input.bam && samtools view -c -F 2304 input.bam > primary_count.txt
```

A positive count and quickcheck do not establish full structural validity. For Picard validation and SN/LN/M5 reference comparison use `references/integrity-dictionary.md`. Compare transfer checksums to detect missing blocks.

Done when the relevant integrity and reference checks are interpreted.
