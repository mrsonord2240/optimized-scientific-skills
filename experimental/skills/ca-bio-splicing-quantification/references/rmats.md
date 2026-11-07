## rMATS-turbo Workflow

**Goal:** Quantify SE/A5SS/A3SS/MXE/RI events from BAMs aligned with STAR 2-pass.

**Approach:** Group BAMs by condition (each `--b1`/`--b2` file is one line of comma-separated BAM paths, e.g. `a.bam,b.bam`), run rMATS with `--statoff` for quantification only, then parse JC.txt files for per-replicate PSI.

```bash
rmats.py \
    --b1 condition1_bams.txt \
    --b2 condition2_bams.txt \
    --gtf annotation.gtf \
    -t paired \
    --readLength 150 \
    --variable-read-length \
    --libType fr-firststrand \
    --nthread 8 \
    --od rmats_output \
    --tmp rmats_tmp \
    --novelSS \
    --statoff
```

Key flags: `--novelSS` discovers junctions absent from the GTF (recommended with STAR 2-pass output). `--variable-read-length` allows mixed read lengths in the cohort. `--libType fr-firststrand` matches Illumina TruSeq stranded; verify with RSeQC `infer_experiment.py`. `--statoff` is for quantification-only runs; omit for differential testing.

```python
import sys
sys.path.insert(0, 'scripts')   # this Skill's scripts/ directory
from quantify_splicing import parse_rmats_output

# SE, A5SS, A3SS, MXE or RI; JC (default) or JCEC via counts='JCEC'.
# mean_PSI pools IncLevel1 + IncLevel2 only (never IncLevelDifference); an event is kept
# when every replicate in BOTH groups has IJC+SJC >= min_junction_reads.
reliable = parse_rmats_output('rmats_output', 'SE', min_junction_reads=20)
```

IncLevel1/IncLevel2 and the IJC/SJC columns are comma-separated per-replicate strings (`NA` where IJC+SJC = 0), so they cannot be averaged directly; `IncLevelDifference` is the group contrast (group 1 minus group 2). With `--statoff`, PValue and FDR are `NA`.

### JC vs JCEC files

`SE.MATS.JC.txt` uses **only junction-spanning reads**. `SE.MATS.JCEC.txt` adds **reads contained within the alternative exon body** as inclusion evidence.

- Prefer **JC** for clean cassette-exon analysis when reads spanning junctions are sufficient.
- Use **JCEC** when alternative exons are short (<50nt) and junction-spanning reads are scarce.
- Avoid **JCEC** when intron retention overlaps the alternative exon body — exon-body reads may come from retained introns, not inclusion isoform.
