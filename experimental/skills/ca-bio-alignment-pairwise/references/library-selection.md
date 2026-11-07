# Pairwise Library Selection (reference)

## Pairwise Library Selection

`Bio.Align.PairwiseAligner` is the right default for interactive use, scripting, and pair sizes up to a few thousand residues, but it is not the fastest or most scalable option. For high-throughput screens, very long sequences, or production pipelines, switch to a SIMD-accelerated or specialised library.

| Library | Speed vs Bio.Align | Alphabet | Scoring | Vectorization | When to use |
|---------|-------------------|----------|---------|---------------|-------------|
| `Bio.Align.PairwiseAligner` (BioPython) | 1x baseline | DNA / RNA / protein | Matrix + affine | C-backed Gotoh | Default, <10 kb pairs, interactive use |
| `parasail` (Daily 2016 BMC Bioinf) | 2-10x (300 nt), ~3-8x (20 kb); measured | DNA / protein | Matrix + affine | SSE / AVX SIMD | High-throughput SW or NW; benchmark loops |
| `edlib` (Sosic & Sikic 2017 Bioinf) | 15-26x (300 nt), 400x+ (20 kb); measured | any alphabet | Edit distance only | Bit-parallel Myers | Read mapping, k-mer search, primer placement |
| `pywfa` / WFA2 (Marco-Sola 2021 Bioinformatics 37:456; BiWFA: Marco-Sola 2023 Bioinformatics 39:btad074) | Best for low-divergence | DNA | Matrix + affine | Wavefront, O(s) memory | Long, near-identical sequences (>10 kb, <5% diverged) |
| `mappy` / minimap2 (Li 2018 Bioinf) | Production reads-to-genome | DNA | Chain + base-level | k-mer chain | Long-read mapping, splice-aware DNA |
| `Bio.pairwise2` | DEPRECATED | -- | -- | -- | Migrate to `PairwiseAligner` (deprecated in BioPython 1.80; not yet removed; migrate proactively) |
| EMBOSS `needle` / `water` | ~Bio.Align | DNA / protein | Matrix + affine | None | Reproducibility, audit trails (fixed, documented default parameters) |

Speed numbers for parasail and edlib were measured (Windows, score only, 5% divergence for 300 nt pairs, 3% for 20 kb; timings vary run to run); the rest are literature figures. Benchmark on representative inputs before committing. Critical caveats:
- **WFA / BiWFA**: 10-100x faster than Gotoh below 5% divergence; above ~10% it converges to Gotoh complexity. Right tool for PacBio HiFi self-similarity or assembly-vs-reference; not for distant homologs.
- **edlib**: measured 15x vs Biopython Levenshtein and 26x vs affine scoring on 1000 pairs of 300 nt, 400x+ on a 20 kb pair; the speedup shrinks as divergence grows (literature: ~64x above ~50% divergence). Its scoring is unit-cost edit distance only (no matrix, no affine gaps; `task='path'` still returns the alignment and CIGAR), so for high-divergence DNA (<70% nucleotide identity) prefer parasail's SIMD score-only mode.
- **parasail**: SIMD only realises its advantage on long sequences in amortised batch loops. **Fixed-width variants silently saturate**: `nw_striped_16` on a 20 kb pair returned score 0 with `.saturated == True` (true score 35565). Use the `*_sat` variants (they widen 8 -> 16 -> 32 bit) and check `.saturated`.

Verified snippets (parasail's open/extend use the Biopython convention, so `10, 1` = `open_gap_score=-10, extend_gap_score=-1`):

```python
import parasail, edlib
mat = parasail.matrix_create('ACGT', 2, -1)                  # protein: parasail.blosum62
r = parasail.nw_striped_sat(query, target, 10, 1, mat)       # sw_striped_sat for local
assert not r.saturated
print(r.score)   # == PairwiseAligner(mode='global', match_score=2, mismatch_score=-1, open_gap_score=-10, extend_gap_score=-1).score(target, query)

# edit distance == -score of PairwiseAligner(match 0, mismatch -1, all gaps -1); mode 'HW' = query anywhere in target
d = edlib.align(query, target, mode='NW', task='distance')['editDistance']
```

`pywfa` and `mappy` (tested on Linux; their Windows builds failed):

```python
from pywfa import WavefrontAligner
x, o, e = 4, 6, 2                      # WFA2 penalties: a gap of length L costs o + e*L
r = WavefrontAligner(target, span='end-to-end', mismatch=x, gap_opening=o, gap_extension=e)(query)
# r.score == PairwiseAligner(mode='global', match_score=0, mismatch_score=-x, open_gap_score=-(o+e), extend_gap_score=-e).score(target, query)

import mappy
for h in mappy.Aligner('ref.fa', preset='map-ont').map(read):   # presets: sr, map-hifi, asm5, ...
    print(h.ctg, h.r_st, h.r_en, h.mapq, h.cigar_str)
```

When uncertain which algorithm Biopython's aligner selected internally, inspect `aligner.algorithm` after configuration -- it returns the resolved variant ("Needleman-Wunsch", "Smith-Waterman", "Gotoh global alignment algorithm", "Gotoh local alignment algorithm", "Waterman-Smith-Beyer global alignment algorithm", "Waterman-Smith-Beyer local alignment algorithm") for deterministic auditing.

## References

- Sosic M, Sikic M. 2017. Edlib: a C/C++ library for fast, exact sequence alignment. Bioinf 33:1394-1395.
- Daily J. 2016. Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments. BMC Bioinf 17:81.
- Marco-Sola S et al. 2021. Fast gap-affine pairwise alignment using the wavefront algorithm. Bioinf 37:456-463.
- Marco-Sola S et al. 2023. Optimal gap-affine alignment in O(s) space. Bioinformatics 39(2):btad074.
