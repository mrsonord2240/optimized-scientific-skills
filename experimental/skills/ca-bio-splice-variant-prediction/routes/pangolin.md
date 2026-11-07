# Tissue-specific prediction

```bash
python scripts/pangolin_tissue.py chr17 7674292 T C --strand - --fasta genome.fa > tissue_scores.txt
```

Replace coordinates, gene strand and FASTA with the input variant. Four model tissues are heart, liver, brain and testis. For VCF scoring split multiallelic records and use Pangolin -m False; see `references/pangolin.md`. Check REF and use uppercase reference sequence.

Done when requested tissue scores and missing inputs are recorded.
