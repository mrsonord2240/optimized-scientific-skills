# Profiles



**Goal:** Add new sequences to an existing MSA without realigning the entire dataset.

**Approach:** Use MAFFT's `--add` for full-length sequences, `--addfragments` for partial / surveillance / metagenomic reads. The two flags are NOT interchangeable.

| Flag | Use when | Behaviour |
|------|----------|-----------|
| `--add` | New sequences are full-length homologues | Each new sequence is profile-aligned end-to-end |
| `--addfragments` | New sequences are partial (reads, contigs, surveillance amplicons) | Free terminal gaps; new sequences may align to a sub-region of the profile |
| `--addprofile` | Adding an entire pre-aligned profile | Profile-profile alignment |
| `--keeplength` | Output must have the same column count as input MSA | Insertions in new sequences are dropped (key for HMMER/Pfam-style use) |

```bash
# Full-length additions; may extend alignment with new columns
mafft --add new_seqs.fasta existing_alignment.fasta > updated.fasta

# Partial reads or fragments; common for SARS-CoV-2 surveillance
mafft --addfragments new_reads.fasta --keeplength reference_msa.fasta > updated.fasta

# Profile-profile merge
mafft --addprofile second_msa.fasta first_msa.fasta > merged.fasta
```

For HMM-curated families (Pfam, Rfam), `hmmalign --trim --outformat afa profile.hmm new_seqs.fa > out.fasta` is the canonical alternative; it constrains insertions to lowercase columns rather than introducing new alignment columns and is the format expected by downstream HMMER/HHsuite tooling.
