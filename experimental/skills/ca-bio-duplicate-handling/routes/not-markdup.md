# Assays where samtools markdup is wrong

Stop. Do not run `scripts/markdup.sh` (it refuses these anyway). Use the alternative.

| Assay | Why not | Do this instead |
|-------|---------|-----------------|
| Bulk RNA-seq, no UMIs | Duplicates are biology at highly expressed loci | Mark nothing |
| Bulk RNA-seq with UMIs, scRNA (10x, STARsolo, drop-seq) | Same | `routes/umi.md` |
| Amplicon or hotspot panel, no UMI | Every read is a coordinate duplicate; markdup erases the data | `samtools ampliconclip` (alignment-amplicon-clipping Skill) |
| Amplicon panel with UMI | Same | `routes/umi.md` (fgbio consensus) |
| Long-read native (ONT, unamplified HiFi) | No PCR step | Mark nothing |
| PacBio HiFi amplicons | The rare long-read case where marking applies | `routes/other-markers.md` (`pbmarkdup`) |
| Microbiome 16S/ITS | Read counts encode community structure | Mark nothing |

Done when you have told the user which row applies and what was done instead.
