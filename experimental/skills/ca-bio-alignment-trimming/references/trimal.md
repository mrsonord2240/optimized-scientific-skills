# Trimal



Pick a trimAl mode by downstream tool. `-gappyout` for HMM profile builds, `-gappyout` or explicit `-strict` for ML tree input (`-strictplus` is documented by trimAl as optimized for Neighbour Joining), `-automated1` (documented as optimized for ML) only when audit-grade reproducibility is not required (heuristic has changed across point releases).

| Mode | Flag | Description |
|------|------|-------------|
| Automated heuristic | `-automated1` | Picks `gappyout`, `strict`, or `strictplus` based on alignment characteristics |
| Gap-only | `-gappyout` | Automatic gap-fraction threshold from gap distribution |
| Strict | `-strict` | Combined gap + similarity criterion |
| Strict-plus | `-strictplus` | Stricter than `-strict`; optimized for Neighbour Joining per trimAl help |
| Manual gap | `-gt 0.5` | Remove columns with > 50% gaps (set fraction explicitly) |
| Manual similarity | `-st 0.5` | Remove columns with similarity below threshold |
| Manual conservation | `-cons 60` | Keep at least 60% of original columns regardless of other criteria |
| Sequence-quality | `-resoverlap 0.8 -seqoverlap 75` | Remove sequences with poor residue/sequence overlap; thresholds are dataset-dependent -- see below |
| HTML report | `-htmlout report.html` | Visual diff of kept/removed columns |
| Column mapping | `-colnumbering` | Output original column indices preserved |

**Goal:** Trim an MSA with trimAl using a heuristic or manual threshold suited to the downstream tool.

**Approach:** Select an automated mode (`-automated1`, `-gappyout`, `-strict`; `-strictplus` for NJ) for typical use, or compose `-gt`, `-st`, `-cons` thresholds for manual control; export an HTML diff or column index mapping when audit-grade reproducibility is required.

```bash
trimal -in input.fasta -out trimmed.fasta -automated1
trimal -in input.fasta -out trimmed.fasta -gappyout
trimal -in input.fasta -out trimmed.fasta -strictplus -htmlout report.html
trimal -in input.fasta -out trimmed.fasta -gt 0.3 -st 0.5 -cons 60 -colnumbering > columns.txt
```

`-gappyout` is the recommended choice for HMM profile building (HMMER `hmmbuild` benefits from aggressive gap removal). For ML tree input use `-gappyout` or `-strict`; `-strictplus` is NJ-oriented and was the only trimAl mode that lost topological accuracy in a 10-replicate simulated check. `-automated1` is the fallback when characterising the dataset is impractical (it chose `-strict` and removed 61% of columns on one 15-sequence protein alignment -- check retention).

**Reproducibility note:** `trimAl -automated1` selects between `gappyout`, `strict`, and `strictplus` via internal heuristics that have changed between releases. For audit-grade reproducibility, do NOT use `-automated1`; specify the underlying mode explicitly and record `trimal --version` in pipeline manifests.

### trimAl Sequence-Overlap Filtering: Check What Was Removed

`-resoverlap`/`-seqoverlap` remove whole SEQUENCES, not columns: a sequence is kept only if at least `-seqoverlap` percent of its non-gap residues fall in columns whose overlap with the rest of the alignment reaches `-resoverlap`. This is meant to catch fragmentary or partial sequences, but `-resoverlap 0.8 -seqoverlap 75` (trimAl's own help-text example) is not a validated default -- on a real alignment it can also drop full-length, legitimate sequences if enough columns are gappy or divergent. Thresholds are dataset-dependent; treat them as a starting point, not a fixed recommendation. In a simulated 15-sequence protein alignment with two genuine partial contigs, `-resoverlap 0.8 -seqoverlap 75` removed 4 of 15 sequences, not 2: the two fragments plus two full-length sequences.

trimAl does not print which sequences it removed -- stdout only lists all-gap columns dropped afterward. Always diff the headers before vs after to see what actually left the alignment:

```bash
trimal -in input.fasta -out filtered.fasta -resoverlap 0.8 -seqoverlap 75
grep '^>' input.fasta | sort > in.txt
grep '^>' filtered.fasta | sort > out.txt
diff in.txt out.txt
```

If a removed sequence is full-length rather than a fragment you meant to drop, the threshold is too strict for this dataset: relax `-seqoverlap` and re-run the diff. On the same 15-sequence alignment, lowering `-seqoverlap` from 75 to 60 dropped only the two genuine fragments and kept both full-length sequences that 75 had discarded.
