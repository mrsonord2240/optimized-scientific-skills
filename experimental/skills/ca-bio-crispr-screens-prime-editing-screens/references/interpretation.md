## Prime Editor Chemistry Comparison

| Editor | Year | Mechanism | Indel rate | Use when |
|--------|------|-----------|------------|----------|
| PE2 (Anzalone 2019) | 2019 | nCas9-RT fusion + pegRNA | 1-3% | Standard PE; lowest indel rate |
| PE3 | 2019 | PE2 + nick of opposite strand by additional sgRNA | 2-5% | Higher editing efficiency, slightly more indels |
| PE3b | 2019 | PE3 with edit-blocking ssgRNA | 1-3% | When PE3's added nick risks unwanted indels |
| PEmax (Chen 2021) | 2021 | Engineered RT + nCas9 | 1-2% | Higher editing rate per pegRNA |
| PE5max (Chen 2021) | 2021 | PE3 plus MMR inhibition (MLH1dn) on the PEmax architecture | 1% | Highest efficiency at favorable sites |
| PE6 / dual-pegRNA (2023) | 2023 | Engineered compact PE; twin-pegRNA systems | Variable | Specific applications |

**Decision rule:** For pooled screens at scale, PE2 or PEmax (single-guide architecture) is preferred over PE3, whose additional nicking sgRNA complicates library architecture. For specific high-efficiency edits, PEmax + PRIDICT2-optimized pegRNA.

## pegRNA Architecture

A pegRNA contains four critical elements that determine efficiency. **The 3' extension is
RTT-then-PBS, not PBS-then-RTT** -- PBS is the pegRNA's 3'-terminal element. Getting this order
backwards is not cosmetic: it makes real CRISPResso2 silently report 0% editing on a true 40%
sample (see the CRISPResso2 section below), because `--prime_editing_pegRNA_extension_seq` never
matches the edited allele under the wrong order.

```
5'  SPACER (20 nt)  -- standard sgRNA spacer; defines target locus via NGG PAM
    +
    SCAFFOLD (~80 nt) -- canonical or recoded scaffold (Chen 2021 recodes it to cut scaffold-incorporation byproducts)
    +
    RTT (Reverse Transcription Template, 10-30 nt) -- encodes intended edit; copied by RT.
        Spans the genomic region DOWNSTREAM of the nick (through the PAM and beyond) --
        this is the window PE actually replaces.
    +
    PBS (Primer Binding Site, 8-15 nt) -- the pegRNA's 3'-terminal element; primes reverse
        transcription off the nicked strand. Complementary to the protospacer sequence
        UPSTREAM of the nick -- it must never include PAM bases.
3'
```

Both RTT and PBS are written on the pegRNA as the **reverse complement** of their genomic
template regions (verified against the real PRIDICT2 CLI's own `RTrevcomp`/`PBSrevcomp`/`pegRNA`
columns and against real CRISPResso2 output -- see below).

**Worked example** (hand-derived, independently confirmed by a real CRISPResso2 2.3.4 run
recovering the planted 40.0% Prime-edited / 33.33% Modified-within-Reference exactly): spacer
`ACGTTGACCTGGAACGTTCA`, NGG PAM immediately 3' of it, C>T edit 11 nt downstream of the nick.
RTT genomic window (18 nt, downstream of nick, edit installed) = `TCATGGCGATCTGTAAGC`; its
reverse complement = `GCTTACAGATCGCCATGA`. PBS genomic window (13 nt, upstream of nick) =
`TGACCTGGAACGT`; its reverse complement = `ACGTTCCAGGTCA`. Extension = RTT-revcomp + PBS-revcomp
= `GCTTACAGATCGCCATGAACGTTCCAGGTCA` -- this is the value that belongs in
`--prime_editing_pegRNA_extension_seq`.

**Key design parameters:**
- **PBS length:** 11-13 nt typical; longer for high-GC contexts; PBS GC fraction critical (35-65% target)
- **RTT length:** 10-20 nt typical; longer for distant edits (10+ bp away from cut)
- **RTT-edit position:** intended edit at position 4-30 from cut site
- **Scaffold:** standard sgRNA scaffold OR the Chen 2021 recoded scaffold, which removes homology with the genomic target and cuts scaffold-derived byproducts. The recoding does not itself raise editing efficiency; that comes from MLH1dn and the PEmax architecture.


## Cas9 vs BE vs PE for Variant Installation

| Approach | Bystander | Indels | Coverage | When to use |
|----------|-----------|--------|----------|-------------|
| Cas9 + HDR | None | High | Variable (depends on template integration) | Precise edits at scale; high indel byproduct |
| Base editor | YES | Low (<5%) | Limited by editing window | C->T or A->G at editable position |
| Prime editor | NONE | Low (<3%) | NGG-PAM within 30 nt of edit | Precise variants; multi-base; transversions |
| Cas9 (no template) | NONE | 70%+ | Anywhere with NGG | LoF only; not variant-specific |

**Decision tree:**
- C->T or A->G at editing-window position: BE (higher efficiency than PE)
- Multi-base / transversion / out-of-window: PE
- LoF without specifying variant: Cas9
- Random insertions: HDR (lower throughput than PE)
- iPSC / primary-cell variant: PE (no bystander confounding)
- Cancer-line variant scanning or drug-resistance variant: PE or BE; cross-validate both (chemistry-dependent) (see `references/pooled-screen-methods.md`)

## Quantitative Thresholds

| Threshold | Value | Source / Rationale |
|-----------|-------|--------------------|
| PRIDICT2 efficiency for library inclusion | >50% | Project-chosen cutoff; PRIDICT2 prescribes none |
| Intended edit % for screen power | >5%; >20% at favorable sites | Field convention |
| Scaffold incorporation | <2% (clean PE); <5% acceptable | Empirical |
| Indel byproduct | <3% (PE2); <5% (PE3) | Anzalone 2019; Chen 2021 |
| PBS GC content | 40-55% | PRIDICT2 |
| PBS length | 11-13 nt | PRIDICT2 |
| RTT length | 10-20 nt | PRIDICT2 |
| Edit position from cut | 1-30 nt | Anzalone 2019 |
| Cell line for PE | K562, HEK293T, U2OS validated | High RT expression |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Low editing across library | Cell-line RT inactivity | Verify PE2 expression; switch to validated line |
| Scaffold incorporation >10% | RTT too short | Re-design with longer RTT |
| Partial multi-base edits | RT processivity limit | Shorter RTT or PE3 |
| PRIDICT predicts but observes much lower | Chromatin context | Score the locus with ePRIDICT (`references/epridict-chromatin.md`); flag low-percentile targets; pilot before library order |
| Library missing variants | No NGG PAM | SpRY-PE; BE alternative |
| CRISPResso2 reports ~0% Prime-edited on a library that should edit | pegRNA extension built PBS-then-RTT instead of RTT-then-PBS | Rebuild extension as RTT-revcomp + PBS-revcomp |
| PRIDICT2 batch: `--summarize: expected one argument` | Bare `--summarize` flag | Pass a value: `--summarize K562` (or `HEK`) |
| PRIDICT2 batch: summary file is `""` (exit 0) | Wrong header (`sequence` not `editseq`), zero data rows, or no NGG PAM near any edit | Read the run's stdout for the cause; see Failure Modes, "Batch run exits 0 with an empty summary file" |
| PRIDICT2 batch: `FileNotFoundError` referencing `input/<file>` | CSV not placed in the CLI's default `./input` directory | `mkdir input` and move the CSV there, or pass `--input-dir` explicitly |
| PRIDICT2 batch `--summarize`: `FileNotFoundError` on the output directory itself | `--output-dir` does not exist yet -- `--summarize` lists existing `.csv` files there before the run starts | `mkdir -p <output-dir>` before running with `--summarize` |
| PRIDICT2 batch `--summarize`: `ValueError: Output directory is not empty. Please move or delete existing .csv files` | `--output-dir` exists but already holds a `.csv` -- including a previous PRIDICT2 run's own per-sequence `_pegRNA_Pridict_full.csv` files, which the run itself writes there | Pre-creating the directory is necessary but not sufficient: a second `--summarize` run into the same `--output-dir` always fails. Point it at a fresh directory (or clear the `.csv` files first) |
| PRIDICT2 batch with `--cores` above 3 | The tool documents 3 as the maximum ("Maximum 3 cores to prevent memory issues") but passes `--cores` straight through unvalidated, so a larger number is accepted silently | Use `--cores 3` (also the default); do not raise it |
| PE concordant with BE on transitions, disagrees on transversions | PE handles transversions BE doesn't | Expected; trust PE |
