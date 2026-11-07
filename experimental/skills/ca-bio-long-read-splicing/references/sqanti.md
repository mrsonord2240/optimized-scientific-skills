## SQANTI3 Classification

**Goal:** Classify discovered isoforms relative to reference; flag artifacts (intra-priming, RT-switching).

**Approach:** Run `sqanti3_qc.py` on the isoform GTF; review classification (FSM/ISM/NIC/NNC/antisense/genic/intergenic/fusion) and quality flags.

```bash
# isoforms.gtf: transcript and exon rows only (FLAIR isoforms.gtf as is; from IsoQuant:
#   awk -F'\t' '$3=="transcript" || $3=="exon"' <prefix>.transcript_models.gtf > isoforms.gtf)
# -o is a file prefix, -d the output directory. --CAGE_peak and --polyA_motif_list are optional: drop them without files
sqanti3_qc.py \
    --isoforms isoforms.gtf \
    --refGTF gencode.v45.annotation.gtf \
    --refFasta reference.fa \
    -o sqanti3 -d sqanti3_qc \
    --aligner_choice minimap2 \
    --CAGE_peak hg38.cage_peak_phase1and2combined_coord.bed \
    --polyA_motif_list human.polyA.list.txt \
    --cpus 8
# classification: sqanti3_qc/sqanti3_classification.txt

sqanti3_filter.py rules \
    --sqanti_class sqanti3_qc/sqanti3_classification.txt \
    --filter_gtf isoforms.gtf \
    -o sqanti3_filtered -d sqanti3_filtered \
    --skip_report
# filtered GTF: sqanti3_filtered/sqanti3_filtered.filtered.gtf; per-isoform filter_result: *_RulesFilter_classification.txt
```

Optional inputs: CAGE peaks (hg38) from refTSS (reftss.riken.jp) or the FANTOM5 file the SQANTI3 docs link (`Magdoll/images_public`, `SQANTI2_support_data/hg38.cage_peak_phase1and2combined_coord.bed.gz`, gunzip first); the polyA motif list from the same repo (`SQANTI2_support_data/human.polyA.list.txt`) or `data/polyA_motifs/` of the SQANTI3 GitHub repo (the conda package ships neither). Both ran on real hg38 isoforms (`within_CAGE_peak` and `polyA_motif_found` filled). `--report skip` on the QC step avoids the HTML/PDF report (the default report also ran, 33 s). `--skip_report` on the filter step is needed on small inputs: its R report failed (tidyselect "subscript out of bounds", exit 1) on 6 and 10 isoforms after the filter tables were written. SQANTI3 6.x has no `--skipORF`.

| SQANTI category | Meaning |
|-----------------|---------|
| FSM (Full Splice Match) | All junctions match reference |
| ISM (Incomplete Splice Match) | Subset of reference junctions |
| NIC (Novel In Catalog) | Novel combination of known junctions |
| NNC (Novel Not in Catalog) | Contains novel junction |
| Antisense | Overlaps gene on opposite strand |
| Genic | Within gene but no junction match |
| Genic intron | Entirely within an intron of a reference gene |
| Intergenic | Between genes |
| Fusion | Spans multiple genes |

SQANTI3 (Pardo-Palacios 2024 *Nat Methods* 21:793-797) is the long-read isoform-curation/QC tool, with structural categories and QC tailored to ONT/PacBio error patterns. **Filter intra-priming and RT-switching** flags before reporting.
