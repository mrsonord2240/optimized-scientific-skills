## Base Editor Chemistry Selection

| Editor | Reaction | Editing window | Indel byproduct rate | When to use |
|--------|----------|----------------|----------------------|-------------|
| BE3 (Komor 2016) | C->T (also G->A on opposite strand) | Pos 4-8 from PAM-distal end | 5-10% | Original; superseded |
| BE4 / BE4max (Koblan 2018) | C->T | Pos 4-8 | <5% | CBE standard |
| eA3A-BE3 | C->T narrow specificity | Pos 5-7 | <5% | Specifically TC contexts (eA3A prefers TC) |
| ABE7.10 (Gaudelli 2017) | A->G (T->C opposite strand) | Pos 4-7 | <2% | First ABE; slow at non-TA contexts |
| ABE8.20 (Gaudelli 2020) | A->G | Pos 4-8 | <2% | Modern ABE; high activity |
| ABE8e (Richter 2020) | A->G | Pos 4-8 | <2% | Highest editing activity; more processive than ABE7.10 |
| evoCDA-BE | C->T (broader) | Pos 1-9 | 5-10% | Larger editing window; more bystander |
| CGBE1 (Kurt 2021) | C->G | Pos 5-7 | 5-10% | C-to-G transversion; rare use |
| GBE (Zhao 2021) | C->G or C->A | Pos 4-7 | 5-10% | Transversions; less mature |

**Decision rule:** For a target SNV at position 4-8 of a candidate spacer with no bystander Cs/As in the same window, BE3-BE4 or ABE7.10 is sufficient. For high-throughput variant scanning where bystander tolerance must be minimized, use eA3A-BE3 (TC contexts only) for C->T, or ABE7.10 rather than ABE8e/ABE8.20 for A->G -- its 4-7 window is the narrowest ABE.

## Editing Window Math

**Why this matters for postdoc-level use:** Base editors are tethered to dCas9 (or nCas9) and the deaminase acts on the displaced ssDNA "R-loop" formed when Cas9 binds. The deaminase has a fixed reach -- positions 4-8 from the PAM-distal end of the protospacer for canonical BE3/BE4, and 4-7 for ABE7.10. Outside this window, editing efficiency drops by 10-50x.

```
PAM-distal end                                                            PAM-proximal
   |                                                                          |
   1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20    NGG
                  ^^^^^^^^^^^
                  Canonical editing window (positions 4-8)

   For BE4max: positions 4-8 are 5-50x more efficient than positions 1-3 or 9-13 (ABE7.10: 4-7)
   For SpABE8e: positions 4-8 (Richter 2020), matching the corresponding CBEs rather than ABE7.10's narrower 4-7
   For evoCDA-BE: window 1-9 (broader; more bystander)
```

**Critical implication for variant interpretation:** If the intended edit is at position 5 and there is an additional editable C/A at position 7, both will be edited in the same molecule. The screen scores the *combination* of edits, not the intended one alone. This is bystander confounding.

## Cas9 vs Base Editor vs Prime Editor for Variant Installation

| Approach | What it does | Bystander | Indels | When to use |
|----------|--------------|-----------|--------|-------------|
| Cas9 + HDR template | Installs precise edit + template | None | High (NHEJ competition) | When precise edit needed; high indel byproduct |
| Cas9 (no template) | Random indels at cut site | None | 70%+ | Loss-of-function; not variant-specific |
| CBE (BE3/BE4) | C->T at editing window | Yes (multiple Cs) | <5% | C->T variants with manageable bystanders |
| ABE (ABE7.10/ABE8e) | A->G at editing window | Yes (multiple As) | <2% | A->G variants; clean for single-A spacers |
| CGBE / GBE | C->G or C->A | Yes | 5-10% | Transversions; rare use cases |
| Prime editor (PE2/PE3) | Templated edit; any base change | None | 1-3% | Precise variants; lower efficiency |

**Decision:** For C->T or A->G with available editing window: base editor is preferred (higher efficiency than PE). For other transitions/transversions, multi-base edits, or zero-bystander requirements: prime editor.
