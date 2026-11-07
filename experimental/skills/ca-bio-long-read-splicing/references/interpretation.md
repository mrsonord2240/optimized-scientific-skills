## Per-Tool Failure Modes

### minimap2: -uf on Unoriented Reads

**Trigger:** `-uf` on ONT cDNA or unoriented HiFi (the recipes above give the rule).

**Symptom:** ~half the reads with false junctions, inferred gene strand wrong for ~half, few reads kept by `flair correct`; the orientation check gives ~0.5.

**Fix:** Drop `-uf` (the presets `splice:hq` and `splice -k14` differ little; see the alignment section).

### IsoQuant: Memory Pressure

**Trigger:** Atlas-scale cohort or low-RAM environment.

**Mechanism:** IsoQuant builds graph structures across all reads simultaneously.

**Symptom:** OOM kill; very slow runtime.

**Fix:** Increase RAM to >=64 GB; or batch by chromosome.

### Bambu: NDR Mistuning

**Trigger:** NDR=0.5+ or NDR=0.01.

**Mechanism:** NDR controls the precision-recall tradeoff for novel transcripts.

**Symptom:** Too many spurious novel transcripts (high NDR) or missing real novel transcripts (low NDR).

**Fix:** Default NDR=0.1 is balanced; adjust based on validation expectations.

### SQANTI3: RT-Switching Flags

**Trigger:** PacBio/ONT cDNA libraries with template switching artifacts.

**Mechanism:** RT-switching produces chimeric reads spanning two unrelated transcripts; SQANTI3 flags these.

**Symptom:** Many "fusion" transcripts in non-cancer samples; biologically implausible.

**Fix:** Filter out RT-switching flags via `sqanti3_filter.py`; investigate library prep if rate >5%.

### FLAIR: Novel Splice Sites Need Orthogonal Junctions

**Trigger:** Running `flair correct` with only `--gtf` (no `--junction_tab`/`--junction_bed`).

**Mechanism:** FLAIR keeps only junctions supported by the annotation or by the orthogonal (short-read) junctions; a novel splice site is moved to the nearest supported one or the read is dropped as inconsistent.

**Symptom:** A planted 24-nt novel-acceptor isoform (20 reads in one sample) was missing from `flair collapse`, and the annotated isoform's count was inflated by exactly those reads (194 vs 174).

**Fix:** Add the short-read junctions (STAR `SJ.out.tab` via `--junction_tab`, or a BED via `--junction_bed`, a 6-column BED: chrom, start, end, name, score = read support, strand). With them all 7 planted isoforms were recovered with exact counts.

## Reconciliation: When Long-Read Tools Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| FLAIR has more isoforms than IsoQuant | FLAIR collapse less stringent; or IsoQuant filtered more aggressively | Both tools have valid pipelines; report based on use case |
| Bambu calls fewer novel than IsoQuant | Bambu NDR=0.1 is more conservative | Adjust NDR or trust Bambu's calibration |
| SQANTI3 classifies as NNC, FLAIR thinks FSM | GENCODE version mismatch | Verify both tools use same annotation |
| Long-read isoform calls don't match short-read events | Short-read EM ambiguity; or long-read coverage gap | Trust long-read for unambiguous; trust short-read for high-coverage events |

## Quality Control for Long-Read Splicing

| Metric | PacBio HiFi | ONT cDNA R10.4.1 |
|--------|-------------|-------------------|
| Read accuracy (modal) | Q30+ (>=99.9%) | ~98% simplex / ~99% duplex |
| Splice junction concordance to short-read truth (unsourced estimate) | ~98% | 95-98% |
| Median read length (transcripts) | 1-4 kb | 0.5-3 kb |
| Throughput per run | ~25M HiFi reads | Tens of millions |
| Library input | 100-500 ng total RNA | 100-500 ng |
| Read direction | TSO + dT primed | TSO or random hexamer |

Figures are approximate vendor or literature values, not measured here.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `Bambu: Input annotation file not readable. Requires .gtf/.gff format or TxDb object` (`prepareAnnotations`) | Empty GTF or rows with fewer than 9 tab-separated columns (reproduced with rows cut to 6 columns; `gffread -E` did not flag them) | `awk -F'\t' '!/^#/ && NF!=9' annotation.gtf` lists the bad rows |

## Quality Thresholds

| Metric | Recommendation | Source |
|--------|----------------|--------|
| Full-length non-chimeric (FLNC) % | >=80% (PacBio Iso-Seq) | PacBio convention |
| FSM% | >=50% in well-annotated genome (field-convention rule of thumb; not specified in the SQANTI paper; real FLAIR-test isoforms gave 7 of 36 FSM, so it depends on annotation and sample) | SQANTI3 documentation; Tardaguila 2018 *Genome Res* 28:396 |
| NNC% | <=30% (>30% suggests artifacts unless biologically interesting) | SQANTI3 convention |
| Junction support | >=2 reads (or >=3 with strict filtering) | Conservative |
| Bambu NDR | 0.1 default; 0.05 stringent | Chen 2023 *Nat Methods* 20:1187 |
| SQANTI3 RT-switching flag | filter out unless validated | SQANTI3 convention |
| SQANTI3 intra-priming flag | filter out | SQANTI3 convention |
| HiFi CCS passes | >=3 | PacBio convention for Q30+ |

## Common Pitfalls

- **Ignoring reference annotation completeness** — SQANTI3 NNC categorization differs by GENCODE version; report version with results.
- **Not running `isoseq refine`** — concatemers and polyA artifacts inflate isoform counts.
- **Confusing FLAIR's 'collapse' with 'cluster'** — collapse merges similar isoforms post-alignment; `isoseq cluster2` merges raw reads pre-alignment.
- **Skipping CAGE / polyA validation in SQANTI3** — TSS / TTS hallucination is common in long-read isoforms.
- **DTU on too few replicates** — long-read is expensive; n=2 vs n=2 is common but underpowered.
