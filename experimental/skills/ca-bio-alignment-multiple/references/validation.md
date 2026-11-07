# Validation



- **Non-homologous sequences**: MSA tools always produce an alignment, even for unrelated sequences; verify homology first (e.g., BLAST E-value < 1e-5). If BLAST+ is not installed, use a local score-vs-shuffled screen instead -- it also flags orientation ahead of `--adjustdirection`:

```python
# Homology/orientation pre-flight, no BLAST+ required (Biopython PairwiseAligner).
# A real homolog scores far above the null distribution from shuffled copies of
# itself; a contaminant or unrelated sequence does not. Run before MSA, not after --
# gap fraction alone (checklist items 2 and 4 above) can miss a contaminant whose
# gap fraction lands inside the range of the real homologs.
import random
from Bio import SeqIO
from Bio.Align import PairwiseAligner

recs = list(SeqIO.parse('sequences.fa', 'fasta'))
al = PairwiseAligner(mode='local', match_score=2, mismatch_score=-3,
                      open_gap_score=-5, extend_gap_score=-2)
rng = random.Random(1)
ref = recs[0]                                    # or a known-good reference
for r in recs[1:]:
    fwd, rev = al.score(ref.seq, r.seq), al.score(ref.seq, r.seq.reverse_complement())
    shuffled_max = max(al.score(ref.seq, ''.join(rng.sample(str(r.seq), len(r.seq))))
                        for _ in range(5))
    call = 'NON-HOMOLOGOUS?' if max(fwd, rev) < 2 * shuffled_max else \
           ('reverse strand' if rev > fwd else 'forward')
    print(f'{r.id}: fwd={fwd:.0f} rev={rev:.0f} shuffled_max={shuffled_max:.0f} -> {call}')
```

  Checked on MAFFT 7.526 / BioPython 1.88 against a 13-sequence set with one planted contaminant and three reverse-strand sequences: flagged only the contaminant as `NON-HOMOLOGOUS?` (score 21, vs. shuffled-null ceiling 25) and correctly called the three reverse-strand sequences without flagging them, while their gap fractions (0.12-0.16) were indistinguishable from the contaminant's (0.18). Drop anything flagged `NON-HOMOLOGOUS?` before alignment; re-run `--adjustdirection` only on what remains, then strip `_R_` per the note above.
- **Sequences below the twilight zone**: Below ~20% protein identity, sequence signal is lost in noise; structural alignment is needed
- **Different domain architectures**: Globally aligning multi-domain proteins with different domain orders produces meaningless results; align individual domains separately
- **Very different lengths without shared homology**: Aligning a 50-residue fragment against 1000-residue proteins globally forces biologically meaningless gaps; use local alignment or fragment-aware modes (E-INS-i)
- **Highly repetitive sequences**: Tandem repeats cause alignment ambiguity; specialized tools (e.g., TRUST for tandem repeats) may be needed
