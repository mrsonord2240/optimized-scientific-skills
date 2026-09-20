# BLAST statistics: E-value derivation, CBS modes, twilight zone

Deeper quantitative background for the two statistical decisions in `SKILL.md`'s "E-value
interpretation" and "Composition-Based Statistics (CBS)" sections. Load this only when you need
the derivation, the full table, or a boundary case -- the two decisions that matter day-to-day
(never compare E-values across databases; use CBS=3 for peptides under 30 aa) are already in
`SKILL.md`.

## E-value interpretation (Karlin-Altschul)

E-value = K * m * n * exp(-lambda * S), where m = effective query length, n = effective database
size, lambda and K are scoring-matrix-dependent constants (Karlin & Altschul 1990 PNAS 87:2264).

| E-value | Bit-score (BLOSUM62, protein) | Interpretation |
|---|---|---|
| < 1e-50 | > 200 | Strong; almost certainly homologous |
| 1e-50 to 1e-10 | 100-200 | Significant; likely homolog |
| 1e-10 to 1e-3 | 50-100 | Marginal; check identity + coverage |
| 0.01 to 10 | 30-50 | Possible remote homolog; needs profile method |
| > 10 | < 30 | Random; not meaningful |

**Key implication of E = K * m * n * exp(-lambda * S):** the same alignment against a 100x larger
database has a 100x larger E-value. Cross-database E-value comparison is meaningless. Bit-score is
database-size normalized and is the right cross-database metric.

For protein remote homology where E is marginal (10^-3 to 10^-1), reach for profile methods:
PSI-BLAST, jackhmmer, HHblits, or Foldseek -- see `remote-homology` skill.

### The "twilight zone" of homology

Below roughly 20-35% pairwise identity, sequence-only search (including BLAST) becomes
statistically unreliable at distinguishing true homology from chance similarity -- the "twilight
zone" (Rost 1999 *Protein Eng* 12:85-94). Below ~20% identity, sequence-only methods should not be
trusted at all; reach for structure-based search (Foldseek) or profile/HMM methods instead -- see
`remote-homology`.

## Composition-Based Statistics (CBS)

Compositional bias inflates significance for low-complexity proteins. The CBS modes (Yu et al.
2006 *Nucleic Acids Res* 34:5966):

| `composition_based_statistics` | Mode | Use when |
|---|---|---|
| 0 | Off | Almost never |
| 1 | F&S 2002 score adjustment | Legacy compatibility |
| 2 | Yu&Altschul 2005 conditional score adjustment | **Default since BLAST+ 2.2.17** -- correct for most cases |
| 3 | Universal statistics | Short queries (< 30 aa) where mode 2 over-corrects |

For protein queries under 30 aa, switch to CBS=3. For protein with known compositional bias (e.g.
coiled-coil regions, signal peptides), CBS=2 is appropriate but consider hard-masking with SEG.
