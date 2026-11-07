# Diagnosing a failed metric

| Symptom | Likely cause | Action |
|---------|--------------|--------|
| chrM over 50% of mapped reads | No nuclear isolation or too little detergent; Tn5 cuts naked mtDNA | Re-prep with Omni-ATAC or fast-ATAC; chrM stripping hides the wasted sequencing |
| NRF under 0.7, PBC2 under 1.0, heavy dedup loss | Over-amplified library, low input (under 5000 cells) with over 12 PCR cycles | No post-hoc fix; re-prep with more cells and fewer cycles. Confirm with PBC1 < 0.7 and pile-ups |
| TSS enrichment under 5 | Over-transposition; genome build mismatch between TSS BED and BAM; flipped strand | Check build and strand column first; profile genuinely flat means over-transposed, lower Tn5 time or amount next prep |
| TSS score 21 against ENCODE 7 | Calculator was ATACseqQC | Recompute ENCODE-style (`routes/tss-encode.md`); the equivalent may be about 8 |
| TSS enrichment low on an Omni-ATAC sample | FFPE-derived TSSs degraded | Use a peak-based metric |
| FRiP under 0.2 with normal peak count | Diffuse signal, or peak set too restrictive | Check TSS enrichment first; if TSS is fine, the caller undercalls: try a looser cutoff (`-p 0.01`) and recompute |
| FRiP differs 2x between identical runs | Peak set drift (q-value cutoff) | Pin caller version and cutoff |
| Insert-size distribution flat after Picard | Single-end sample | Insert size is paired-end only |

Fragment-size patterns (read the PDF from `routes/fragment-tss-frip.md`):

| Pattern | Interpretation | Action |
|---------|----------------|--------|
| NFR (about 50 bp) >> mono (about 200) > di (about 400) > tri (about 600) | Excellent | Pass |
| NFR + mono only, di and tri faint | Acceptable, common in Omni-ATAC | Pass |
| Single broad peak, flat after NFR or no NFR | Over-transposition or degraded chromatin | Reject |
| Mono >> NFR | Under-transposition or condensed chromatin | Caution; low peak counts |
| Sharp 147 bp spike with no flanks | MNase-like contamination | Reject, not ATAC-grade |
| Sub-peaks every 10.4 bp on the NFR | Helical phasing visible | Pass; positive indicator, absence is not a flag |

Done when the cause is named and the keep, re-prep or drop decision is stated.
