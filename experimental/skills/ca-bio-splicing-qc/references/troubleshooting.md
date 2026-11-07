## Per-Tool Failure Modes

### RSeQC junction tools: silent zeros

**Trigger:** BED/BAM contig names differ, the gene model is BED6, or there are no spliced reads at MAPQ >= 30 (`-q`, default 30).

**Symptom:** exit 0 with all-zero saturation curves, `Total 0 usable reads were sampled` (`infer_experiment.py`), an empty `.junction.xls` (0 bytes, `pandas.errors.EmptyDataError`), or 0% known because every junction is `complete_novel` (BED6).

**Fix:** run the checks in **Before You Run Anything**; the helper raises a clear error for each. For `infer_experiment.py` on a BAM whose MAPQ never reaches 30 (multi-mapper-heavy aligner output), pass `-q 0`; `-q 30` is already the default.

### `infer_experiment.py`: "0 usable reads"

`Total 0 usable reads were sampled` / `Unknown data type: Mixture` means the gene model shares no contig with the BAM, the file is a GTF rather than BED, or every read fails `-q`. It is not a sample-size problem; `-s` (default 200000) only caps the sample.

### MaxEntScan: invalid sequences

See Splice Site Strength: wrong length ends the process, N/U raise `KeyError`; validate first.

## Common Errors

| Error (as printed) | Cause | Solution |
|--------------------|-------|----------|
| `EXITING because of fatal PARAMETERS error: present --sjdbOverhang=74 is not equal to the value at the genome generation step =100` | Index built with a different overhang | Use the index's value, or regenerate the index with `--sjdbOverhang` = read length - 1 |
| `Fatal LIMIT error: the number of junctions to be inserted on the fly =10552 is larger than the limitSjdbInsertNsj=...` / `SOLUTION: re-run with at least --limitSjdbInsertNsj 10552` | Annotation + merged novel junctions exceed the default 1,000,000, or a small value was set | Re-run with the value STAR prints; filter the merged file if it is unexpectedly large |
| `junction_annotation.py: error: Rscript executable not found: Rscript` (exit 1) | R not on PATH | `--skip-plot`, or install R |
| `pandas.errors.EmptyDataError: No columns to parse from file` on `.junction.xls` | BAM has no spliced reads at MAPQ >= 30 | Check the library; lower `-q` only if the aligner's MAPQ scale requires it |
| `ValueError: fetch called on bamfile without index` | pysam `fetch()` on an unindexed BAM | `samtools index`, or use `scripts/splicing_qc.py` (reads with `until_eof`) |
| `SystemExit: Wrong length of fa!` / `KeyError` from `maxentpy` | 5'ss not 9 nt, 3'ss not 23 nt, or N/U present | Use `score_splice_sites` from `scripts/splicing_qc.py` |

## Quality Thresholds

| Metric | Good | Acceptable | Poor | Source |
|--------|------|------------|------|--------|
| Read length (PE) | >=150 nt | 100-149 nt | <100 nt | convention |
| Sequencing depth | >=50M | 30-49M | <30M | convention |
| Junction saturation (known curve, growth 80->100%) | <2% (plateau) | Near plateau | Still rising | RSeQC convention |
| Known-junction fraction (reads) | >=80% | 60-80% | <60% (suspect or interesting) | RSeQC convention |
| Junctions with >=10 anchored reads | >=50% | 30-50% | <30% | convention; depends on depth (11% on a 50k-pair chrX sample) |
| Strandedness (dominant direction) | >=90% | 70-90% | <70% (treat as unstranded) | RSeQC convention |

Splice-site MaxEnt, rRNA and Picard cut-offs are in their own sections.

## Troubleshooting Low Event Detection

| Issue | Possible causes | Solutions |
|-------|-----------------|-----------|
| Few events called | Low depth; short reads; SE; wrong strand | Increase depth; use PE150; verify `--libType` |
| High novel junctions | Annotation gaps; mapping artifacts; biology (TDP-43, SF3B1) | Update annotation; check 2-pass; consider biology |
| Low IR detection | poly(A) library | Use rRNA depletion |
| Many weak splice sites | Cryptic splicing | Validate with MaxEnt + SpliceAI |
| PSI variance high across replicates | Library prep / RIN inconsistency | Check RIN; check 3' bias (Picard) |
