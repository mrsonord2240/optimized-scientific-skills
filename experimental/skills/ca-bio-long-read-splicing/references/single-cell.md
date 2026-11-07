## Single-Cell Long-Read for Splicing

**Goal:** Combine cell typing (10X 5' short read) with full-length isoform structure (Kinnex / MAS-Iso-seq).

**Approach:** Split 10X library; sequence half short-read for cell typing, half PacBio Kinnex for isoforms; de-array the Kinnex reads with skera, then run the Iso-Seq steps (lima, isoseq refine, isoseq cluster2). Barcode assignment to cells is outside this Skill (see `single-cell-splicing`).

```bash
# De-array Kinnex/MAS HiFi reads into segmented reads (S-reads). Obtain the
# kit-specific adapter FASTA from PacBio for the exact kit; do not guess its name.
skera split \
    raw_kinnex.bam \
    kit_adapters.fasta \
    segmented.bam

# Then lima -> isoseq refine -> isoseq cluster2 (bioconda lima, isoseq; see long-read-sequencing/isoseq-analysis)
```

skera 1.4.0 was run on a synthetic array. With a valid `zm` tag on its toy BAM, `skera -> lima --isoseq -> isoseq refine` completed (132 input S-reads -> 132 lima -> 132 FLNC). Without `zm`, the toy read names used `?` and lima stalled; this was a fixture defect, not a demonstrated lima limitation. No real Kinnex data was run.

Joglekar et al 2024 (*Nat Neurosci* 27:1051-1063) used this approach to map single-cell isoforms across developing and adult mouse and human brain. See `single-cell-splicing` for tools that work on the demultiplexed data.
