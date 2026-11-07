# Which method

| Scenario | Use | Why |
|----------|-----|-----|
| 16S V4 (about 253 bp), 2x250 | DADA2 paired (`routes/dada2-16s.md`) | Large merge slack; truncate to quality |
| 16S V3-V4 (about 460 bp), 2x250 | DADA2 paired, protect the overlap, loosen reverse `maxEE` | About 28 bp slack |
| ITS | cutadapt, ITSxpress, DADA2 with `truncLen 0` (`routes/its.md`) | Fixed truncation slices real biology |
| Full-length 16S, PacBio CCS | `qiime dada2 denoise-ccs` (`routes/qiime2.md`) | Single-end, resolves to species or strain |
| Several sequencing runs | Per-run inference, `mergeSequenceTables`, one chimera removal | Error model is per run |
| Fast, fixed length, many runs, 16S only | Deblur (`routes/qiime2.md`) | Per-sample independent |
| Rare or singleton ASV sensitivity | DADA2 `--pool pseudo` | Approximates full pooling in linear time |
| NovaSeq, NextSeq, iSeq (binned quality) | Inspect `error_fit_*.png`, enforce a monotonic fit | Few quality bins starve the loess fit |
| 97% OTUs required (legacy) | VSEARCH clustering | The OTU path |
| Shotgun reads | `routes/other-skills.md` | Different category |

ASV versus OTU: ASVs are exact sequences, portable across studies and the field default. They can over-split a genome with divergent 16S copies, so treat ASV count as an upper bound on richness and collapse to a taxonomic rank before richness claims. Species-level 16S calls are usually overconfident, and reads are not cells (copy number 1-15+).

Done when one method is chosen and its route opened.
