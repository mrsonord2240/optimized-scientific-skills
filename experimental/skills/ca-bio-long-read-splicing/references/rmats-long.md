## rMATS-long for Differential Isoform Analysis on Long-Read Data

**Goal:** Apply differential isoform analysis to long-read transcript abundance with classification and visualization.

**Approach:** rMATS-long is a multi-script Python pipeline distributed via bioconda; entry point is `rmats-long` followed by the script name. It supports two modes: **abundance-based** (using ESPRESSO-style abundance estimates) and **ASM-based** (Alternative Splicing Modules — sets of isoforms sharing exon-junction structure). Run preprocessing scripts in order before `rmats_long.py`.

```bash
# conda install -c conda-forge -c bioconda rmats-long   (own env)
set -euo pipefail   # the steps below otherwise carry on after a failure and write header-only tables
mkdir -p alignment_info

# Preprocessing pipeline (ASM mode); per-script flag names verified vs Xinglab/rmats-long
rmats-long organize_gene_info_by_chr.py --gtf annotation.gtf --out-dir gene_info_by_chr/

# simplify_alignment_info processes one sorted, indexed BAM at a time -> one TSV;
# samples.tsv (sample_id<TAB>tsv_path) is what organize_alignment_info_by_gene_and_chr.py reads
: > samples.tsv
for bam in *.bam; do
    id="${bam%.bam}"
    rmats-long simplify_alignment_info.py --in-file "$bam" --out-tsv "alignment_info/${id}.tsv"
    printf '%s\talignment_info/%s.tsv\n' "$id" "$id" >> samples.tsv
done
# every per-sample table must exist and be non-empty
while IFS=$'\t' read -r id tsv; do [ -s "$tsv" ] || { echo "empty or missing $tsv" >&2; exit 1; }; done < samples.tsv

rmats-long organize_alignment_info_by_gene_and_chr.py \
    --gtf-dir gene_info_by_chr/ \
    --out-dir organized/ \
    --samples-tsv samples.tsv

rmats-long detect_splicing_events.py --align-dir organized/ --gtf-dir gene_info_by_chr/ --out-dir events/
rmats-long create_gtf_from_asm_definitions.py --event-dir events/ --out-gtf asm.gtf
rmats-long count_reads_for_asms.py --align-dir organized/ --event-dir events/ --gtf-dir gene_info_by_chr/ --out-dir asm_counts/

# Main differential analysis (ASM mode)
# --group-1 / --group-2 each take the PATH to a file whose single line is a
# comma-separated list of sample IDs (matching the BAM basenames in --align-dir).
echo 'ctrl1,ctrl2,ctrl3' > group1.txt
echo 'trt1,trt2,trt3' > group2.txt
rmats-long rmats_long.py \
    --group-1 group1.txt \
    --group-2 group2.txt \
    --event-dir events/ \
    --asm-counts-dir asm_counts/ \
    --align-dir organized/ \
    --gtf-dir gene_info_by_chr/ \
    --out-dir rmats_long_output/ \
    --adj-pvalue 0.05 \
    --delta-proportion 0.05 \
    --average-reads-per-group 10

# the result tables must have data rows, not just a header
for t in differential_asms.tsv differential_isoforms.tsv; do
    [ "$(wc -l < "rmats_long_output/$t")" -gt 1 ] || { echo "rmats_long_output/$t has no rows" >&2; exit 1; }
done
```

Alternative, abundance-based mode (when you already have ESPRESSO-style estimates), instead of the ASM steps:

```bash
rmats-long rmats_long.py --abundance abundance.esp --updated-gtf updated.gtf \
    --group-1 group1.txt --group-2 group2.txt --out-dir rmats_long_output/ --no-splice-graph-plot
```

Key flags: `--adj-pvalue` (default 0.05), `--delta-proportion` (default 0.05), `--average-reads-per-group` (default 10), `--no-splice-graph-plot` (skip expensive splice-graph rendering).

rMATS-long is a separate tool from short-read rMATS-turbo. The predecessor `lr2rmats` used long reads only to *augment* the short-read rMATS GTF. The ASM framework treats AS as a set-of-isoforms problem, more natural for long-read data than rMATS-turbo's pre-defined event categories.
