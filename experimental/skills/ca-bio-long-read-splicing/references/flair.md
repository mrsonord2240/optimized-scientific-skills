## FLAIR Workflow (correct -> collapse -> quantify -> diffSplice)

**Goal:** Identify, quantify, and test full-length isoforms from long-read RNA-seq across conditions.

**Approach:** Correct splice junctions against annotation and/or short-read junctions, collapse isoforms, quantify per-sample expression, run diffSplice for differential isoform usage.

```bash
# BAM -> BED12 (or `flair align`), one file per sample
bedtools bamtobed -bed12 -i aligned.bam > aligned.bed

# correct takes no genome; needs --gtf and/or short-read junctions (STAR SJ.out.tab via --junction_tab, or BED via --junction_bed)
flair correct \
    --query aligned.bed \
    --gtf gencode.v45.annotation.gtf \
    --junction_tab SJ.out.tab \
    --output flair_corrected \
    --threads 16

# collapse all samples together: cat the corrected BEDs, pass every reads file
flair collapse \
    --query flair_corrected_all_corrected.bed \
    --reads sample.fastq.gz \
    --genome reference.fa \
    --gtf gencode.v45.annotation.gtf \
    --output flair_collapsed \
    --threads 16

# reads_manifest.tsv: no header, tab-separated: sample_id  condition  batch  /path/reads.fastq.gz
flair quantify \
    --reads_manifest reads_manifest.tsv \
    --isoforms flair_collapsed.isoforms.fa \
    --output flair_quantified \
    --threads 16

flair diffSplice \
    --isoforms flair_collapsed.isoforms.bed \
    --counts_matrix flair_quantified.counts.tsv \
    --out_dir flair_diffsplice \
    --test \
    --threads 16
```

FLAIR (Tang 2020 *Nat Commun*) handles ONT and PacBio with the same workflow. `flair quantify` writes `<output>.counts.tsv` (first column `ids` = `<isoform>_<gene>`, sample columns `<id>_<condition>_<batch>`; `--sample_id_only` gives `ID` and plain sample ids). `flair diffSplice` writes per-event-type (`es`, `alt5`, `alt3`, `ir`) inclusion/exclusion tables; `--out_dir` must not exist yet (`-of` overwrites). `--test` adds a DRIMSeq test per event type and needs an `Rscript` on PATH with DRIMSeq, argparse, and data.table; pip/conda FLAIR provides none of them. An event type may have rows but fail DRIMSeq filtering with `No genes left after filtering!`; FLAIR can still exit 0. For every event type required by the analysis, require its `drimseq_<event>_*.tsv` output rather than treating process exit status as success.
