# Bmge Divvier



Use BMGE for deep prokaryotic phylogenomics (the GToTree default) or whenever entropy-based, matrix-aware column filtering is preferred over gap-only heuristics.

**Goal:** Filter MSA columns by matrix-aware entropy and gap-rate thresholds, particularly for deep prokaryotic phylogenomics.

**Approach:** Invoke BMGE with `-t` (sequence type) and tune the entropy and gap-rate thresholds against the divergence depth of the dataset; recalibrate the entropy threshold if the substitution matrix is changed via `-m`. BMGE 1.12 and 2.0 use different flags: on 2.0, `-h` prints help and exits 0 without writing anything, so always check that the output file exists and is non-empty.

```bash
# BMGE 1.12
java -jar BMGE.jar -i input.fasta -t AA  -h 0.5 -g 0.2 -of trimmed.fasta
java -jar BMGE.jar -i input.fasta -t DNA -m DNAPAM100:2 -h 0.5 -g 0.2 -of trimmed.fasta

# BMGE 2.0 (entropy is -e; types AA|CO|NT; output -o or -of)
java -jar BMGE.jar -i input.fasta -t AA -e 0.5 -g 0.2 -o trimmed.fasta
java -jar BMGE.jar -i input.fasta -t NT -m DNAPAM100:2 -e 0.5 -g 0.2 -o trimmed.fasta

test -s trimmed.fasta || echo "BMGE wrote no alignment: check flags against the installed version"
```

| Meaning | BMGE 1.12 | BMGE 2.0 |
|---------|-----------|----------|
| Entropy threshold (lower = more aggressive) | `-h`, default 0.5 | `-e`, default 0.5 |
| Gap rate threshold | `-g`, default 0.2 | `-g`, default 0.5 |
| Substitution matrix | `-m`, default BLOSUM62 (AA) | `-m`, default BLOSUM30 (AA) / DNAPAM180:2 (NT) |
| Sequence type (required) | `-t AA`, `CODON`, `DNA` | `-t AA`, `CO`, `NT` |
| Minimum block size | `-b`, default 5 | `-b`, default 3 |
| Help | -- | `-h` |

Lower entropy thresholds trim much harder: on a simulated deep 16-taxon protein supermatrix, BMGE 1.12 `-h 0.4 -g 0.2` removed 87% of columns, shrank tree length from 12.5 (true) to 8.4 and lowered support, while `-h 0.5` removed 71%. Check retention against the 20%/40% rule and compare trees before lowering the threshold.

**Matrix-aware threshold:** BMGE computes entropy weighted by the chosen substitution-matrix probabilities. Entropy thresholds are calibrated for the version's default matrix (BLOSUM62 in 1.12, BLOSUM30 in 2.0). Switching to `-m BLOSUM30` (deep phylogenomics) makes the same numeric `-h` more permissive; `-m BLOSUM90` (close orthologs) makes it more aggressive. When using a non-default matrix, recalibrate by running on a known-good benchmark before applying.

## Divvier: Column Splitting Instead of Removal

Use Divvier when phylogenetic signal in indels matters and pure removal would discard block-boundary information; it splits ambiguous columns rather than dropping them.

```bash
# Divvier ships as a compiled binary; invoke as ./divvier or `divvier` if on PATH
./divvier -divvy input.fasta

./divvier -partial -mincol 4 -divvygap input.fasta
```

`-divvy` (default) outputs `input.fasta.divvy.fas` (full divvying); `-partial` only filters individual ambiguous characters. `-mincol N` sets the minimum number of confident characters required to keep a split column (integer count, default 2). `-divvygap` writes gaps instead of asterisks for phylogenetics-tool compatibility. Divvier is particularly useful when input alignments have many short conserved blocks separated by ambiguous regions -- removing the regions discards information about block boundaries that splitting preserves.

Maintenance note: Divvier's last release was 2019 (`simonwhelan/Divvier`); the tool is no longer actively maintained but results remain reproducible with the v2019 binary.
