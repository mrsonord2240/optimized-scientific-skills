# Formats and annotations



Three Python libraries cover the alignment-format space, with overlapping but non-identical support. Pick by what is actually required. Checked on Biopython 1.88 / pyhmmer 0.12.3.

| Format | `Bio.AlignIO` (`MultipleSeqAlignment`) | `Bio.Align` (`Alignment`) | `pyhmmer.easel` (`format=`) | Notes |
|--------|---------------|----------------------|-----------------|-------|
| Aligned FASTA | R/W | R/W | R/W (`afa`) | Most portable; loses annotations |
| Clustal | R/W | R/W | R/W | Clustal conservation marks NOT round-tripped |
| PHYLIP (interleaved/sequential/relaxed) | R/W | R/W (strict 10-char names only: no relaxed variant; long ids raise on read and are silently truncated on write) | R/W (`phylip`, `phylips`) | See PHYLIP pitfalls |
| Stockholm | R/W | R (plain files only; TypeError on real Pfam), W (AttributeError) | R/W | Only format preserving GS/GR/GC annotations (GF header lines are dropped by AlignIO) |
| NEXUS | R/W | R/W | -- | MrBayes / PAUP* input; write needs a molecule type |
| MAF (Multiple Alignment Format) | R/W | R/W | -- | UCSC whole-genome alignments |
| A2M | -- (`'fasta'` parser works on padded A2M only) | R (padded A2M), W (needs `column_annotations['state']`, so only for alignments read from A2M) | R/W | HMMER `hmmalign` writes ragged A2M |
| A3M | -- | -- | -- | HH-suite / ColabFold; convert to A2M first (below) |
| MSF (GCG) | R | R | -- | GCG legacy |
| EMBOSS / Mauve XMFA / FASTA-m10 | R (Mauve also W; FASTA-m10 AlignIO only) | R (Mauve also W; no FASTA-m10) | -- | Mostly one-way |
| PSL / chain / BED / SAM / exonerate / bigMaf / bigPsl / bigBed | -- | R/W (pairwise alignments, not MSAs) | -- | Use Kent tools for manipulation |
| HHR / tabular (BLAST) | -- | R only (search hits, not MSAs) | -- | |

**Formats NOT in BioPython** (use dedicated tools):

| Format | Tool | Why |
|--------|------|-----|
| HAL | progressiveCactus, halTools | HDF5-backed multi-genome alignments at TB scale |
| net | UCSC Kent tools (`chainNet`) | Pairwise genome alignment |
| AXT | BLASTZ / lastz native | Pairwise alignment blocks |
| GFA / rGFA | `vg`, `odgi`, `pggb`, gfatools | Pangenome graph |
| GAF | `vg surject`, `vg call` | Graph alignment format (read-to-graph) |

Use `Bio.AlignIO` for MSA files, especially Stockholm and NEXUS; use `Bio.Align` (`Alignment` objects with `.counts()` and `.substitutions`, see "Alternative: Bio.Align Module I/O") when those features are needed. For multi-gigabyte Stockholm databases such as Pfam-A.full, read [Streaming large Stockholm databases](references/stockholm-streaming.md) before choosing an in-memory parser.




Choosing the output format depends on which downstream tool consumes the alignment:

| Downstream Tool | Required Format | BioPython Format String |
|----------------|-----------------|------------------------|
| RAxML-NG, IQ-TREE | PHYLIP (relaxed) | `'phylip-relaxed'` |
| MrBayes | NEXUS (quote-free ids) | `'nexus'` |
| PAUP* | NEXUS or PHYLIP | `'nexus'` or `'phylip'` |
| HMMER, Infernal | Stockholm | `'stockholm'` |
| Pfam/Rfam databases | Stockholm | `'stockholm'` |
| PAML/codeml | PHYLIP (sequential) | `'phylip-sequential'` |
| Most tools | FASTA | `'fasta'` |

### Annotation Preservation

Not all formats support annotations. Converting between formats can silently discard metadata:

| Format | Sequence Annotations | Column Annotations | Secondary Structure |
|--------|---------------------|-------------------|-------------------|
| Stockholm | Yes (GS/GR lines) | Yes (GC lines) | Yes (SS_cons) |
| NEXUS | Partial (SETS block) | Via CHARSET | No |
| Clustal | No (conservation marks not parsed) | No | No |
| PHYLIP | No | No | No |
| FASTA | No | No | No |

Converting Stockholm to FASTA or PHYLIP discards all annotations, secondary structure markup, and per-residue quality scores. Even Stockholm to Stockholm through `AlignIO` drops the `#=GF` header lines (`ID`, `AC`, `DE`, ...; verified on Pfam PF00042). If annotations matter, keep the original Stockholm file as the master copy.

## Format-Specific Notes

### PHYLIP Format Pitfalls

PHYLIP has two incompatible variants (interleaved vs sequential) and two name-length modes (strict vs relaxed). Mixing them up usually fails loudly (`ValueError` on read), but strict names are truncated to 10 characters.

**Strict PHYLIP** truncates sequence names to exactly 10 characters. Writing distinct names that share a 10-character prefix (e.g., `Homo_sapiens_chr1` and `Homo_sapiens_chr2`) raises `ValueError: Repeated name 'Homo_sapie'` (`phylip-sequential` truncates and raises the same way). The silent case is *reading* a foreign strict file whose names collide: both records are accepted with the same id. Names that stay unique are truncated without warning.

```python
# Strict PHYLIP (10-char names, interleaved) -- only for tools requiring it
alignment = AlignIO.read('file.phy', 'phylip')

# Sequential PHYLIP (10-char names, one sequence at a time) -- PAML/codeml.
# Shorten ids to <=10 unique characters first (NCBI headers such as 'lcl|NM_001...' collide)
alignment = AlignIO.read('file.phy', 'phylip-sequential')

# Relaxed PHYLIP (no name limit) -- RAxML-NG, IQ-TREE (recommended default)
alignment = AlignIO.read('file.phy', 'phylip-relaxed')

# Always prefer phylip-relaxed for writing unless the downstream tool
# specifically requires strict format
AlignIO.write(alignment, 'output.phy', 'phylip-relaxed')
```

**Clustal IDs:** Biopython's Clustal writer truncates identifiers to 30 characters without a collision check. After writing Clustal, re-read it and assert the ids remain unique before using it as an input to another tool.

#### PHYLIP-Relaxed Dialect Mismatches Between Tree Tools

Biopython's `'phylip-relaxed'` writes a single space between name and sequence. RAxML-NG and IQ-TREE accept this (checked with names of 138 characters); PhyML and RAxML-NG reject punctuation in names (table below); PAML's codeml expects sequential format with name-truncation behaviour distinct from interleaved. Common silent failures:

| Symptom | Cause | Fix |
|---------|-------|-----|
| RAxML-NG: `ERROR: Invalid character in sequence N at position P: *` | `*` in a nucleotide alignment (RAxML-NG 1.2.2 and 2.0.3 accept `*` in protein data as undetermined, so a stop codon is not scored) | Replace with `N` or `-` |
| IQ-TREE 3.1.3: `WARNING: Some sequence names are changed` | Foreign file with `:` in a name (IQ-TREE renames it; Biopython's relaxed writer already turns `:` into `|` and drops `(` `,`) | Sanitize names yourself: `re.sub(r'[^A-Za-z0-9_]', '_', record.id)` |
| PhyML: `Character ':' is not permitted in sequence name` (same for `,`) | Punctuation in a name (PhyML accepts `(` `)` but writes them into the Newick tree; 138-character names are kept intact, checked on 3.3.20220408 and 3.3.20260528) | Sanitize names with the same `[A-Za-z0-9_]`-only recipe |
| RAxML-NG: `ERROR: Following taxon name contains invalid characters` | `:` `,` `(` or `[` in a name | Sanitize names with the same `[A-Za-z0-9_]`-only recipe |
| codeml: `Error in sequence data file ... separate the sequence from its name by 2 or more spaces` | Used `phylip-relaxed` (single space, interleaved) instead of `phylip-sequential` | codeml requires sequential with short unique names |

IQ-TREE has no validate-only flag (`--check` is invalid). Before a long run, check the input with a zero-iteration run (checked on IQ-TREE 3.1.3; a wrong sequence length gives `ERROR: Line N: Sequence X has wrong sequence length`), or re-read the file with `AlignIO.read`:

```bash
iqtree3 -s file.phy -n 0 -m LG -redo -pre check   # prints 'Alignment has N sequences with M columns'
```

### MAF Block Coordinate Conventions

Read [MAF coordinates](references/maf-coordinates.md) before mapping MAF rows to a genome: minus-strand `start` values require a coordinate conversion and Biopython reports `strand` as an integer.

### Stockholm Format Annotations

Stockholm format (used by Pfam, Rfam, HMMER) supports four annotation line types:

| Line Prefix | Scope | Description | Example |
|-------------|-------|-------------|---------|
| `#=GF` | File | Alignment-level metadata (ID, accession, description) | `#=GF AC PF00001` |
| `#=GC` | Column | Per-column annotation (1 char per alignment column) | `#=GC SS_cons ..(((...)))..` |
| `#=GS` | Sequence | Per-sequence free text (organism, description) | `#=GS seq1 OS Homo sapiens` |
| `#=GR` | Residue | Per-residue annotation (1 char per residue) | `#=GR seq1 SS ..HHH..EEE..` |

Common GC annotations: `SS_cons` (consensus secondary structure), `RF` (reference coordinates), `seq_cons` (consensus sequence).

**WUSS notation** in RNA `#=GC SS_cons` lines uses nested bracket pairs (`<>`, `()`, `[]`, `{}`) for paired bases and characters like `_`, `-`, `,`, `:`, `.`, `~` for unpaired regions; pseudoknots use upper/lower-case letter pairs (`Aa`, `Bb`). Consult the Infernal user guide for the full character table before writing or parsing custom SS_cons strings.

```python
alignment = AlignIO.read('pfam.sto', 'stockholm')

for record in alignment:
    print(record.id, record.annotations)
    if 'secondary_structure' in record.letter_annotations:
        print(f'  SS: {record.letter_annotations["secondary_structure"]}')

ss_cons = alignment.column_annotations.get('secondary_structure')
```

**Round-trip caveat:** re-reading and re-writing as Stockholm preserves GS/GR/GC, but dropped/added sequences invalidate the per-residue annotations -- regenerate annotations after edits.

**Pfam-style `name/start-end` identifier convention:** Pfam, Rfam, and Dfam Stockholm IDs (e.g. `Q9Y6Y0/45-198`) encode a 1-based inclusive region. `AlignIO` already copies the region into `record.annotations['start']` / `['end']` (and the name into `['accession']`) but leaves the suffix on `record.id`. IQ-TREE 3.1.3 accepts ids with `/`; strip the suffix from `record.id` only for a tool that rejects it, then restore it after.

### A2M / A3M Conventions

Read [A2M and A3M conventions](references/a2m-a3m.md) before treating HMMER, HH-suite, or ColabFold output as a rectangular alignment. It covers ragged `hmmalign` output and the first-record rule for `reformat.pl`.

### Streaming Large Stockholm Databases

For Pfam-A.full, BFD, or another collection too large for memory, use [Streaming Stockholm databases](references/stockholm-streaming.md) instead of an in-memory `AlignIO` recipe.
