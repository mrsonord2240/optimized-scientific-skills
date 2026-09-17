# Selection-Analysis Trimming: TCS and MACSE

Moved out of `SKILL.md` (bio-alignment-trimming) because these apply only when the downstream step is
a selection / dN/dS analysis (PAML codeml, HyPhy) -- consult this file when the Goal-Driven Tool
Selection table or Decision Tree in `SKILL.md` points you here.

## TCS Column Masking for Selection Analysis

TCS (Chang et al 2014 MBE) scores each column on a 0-9 scale by consistency across a pairwise-alignment library. The recommended dN/dS workflow:

**Goal:** Mask unreliable columns before selection analysis to prevent alignment errors from inflating false-positive dN/dS calls.

**Approach:** Build a library-rich consistency score via M-Coffee, evaluate the alignment against that library to obtain per-column TCS scores, then keep columns at or above a chosen confidence threshold using `seq_reformat`.

```bash
# 1. Generate library-rich consistency scores via M-Coffee (combines libraries from MAFFT, MUSCLE, ClustalW, ProbCons)
t_coffee input.fasta -mode mcoffee -output fasta_aln -outfile aligned.fasta

# 2. Compute per-column TCS scores; T-Coffee writes aligned.score_ascii itself (not to stdout)
t_coffee -infile aligned.fasta -mode evaluate -output score_ascii

# 3. Filter COLUMNS at a chosen threshold (5-9 = retain columns scoring 5 or higher).
#    +use_cons +keep acts on column scores; +keep alone keeps individual residues within the range.
t_coffee -other_pg seq_reformat -in aligned.fasta -struc_in aligned.score_ascii \
    -struc_in_f number_aln -action +use_cons +keep '[5-9]' -output fasta_aln > aligned_tcs5.fasta
```

Threshold guidance: TCS >= 5 retains columns confidently aligned by the majority of library methods; TCS >= 7 is a stricter operational convention on the seq_reformat 0-9 display scale. Chang et al 2014 do not prescribe an integer 5/7 cutoff -- on BAliBASE 3 and PREFAB 4 structural benchmarks they keep residues scoring above ~0.6 on the native 0-1 scale and drop columns scoring below 2, with TCS outperforming GUIDANCE and HoT. TCS-from-mcoffee gives substantially better column-confidence ranking than TCS-from-default-tcoffee because the library is more diverse. For the PAML branch-site test, Fletcher & Yang 2010 showed alignment errors inflate false positives dramatically (up to ~0.99 with ClustalW, ~0.13 even with codon-aware PRANK) and that removing gappy columns did not reduce them; confidence-based column masking (TCS, GUIDANCE2) is a common mitigation but was not evaluated in that study, since TCS postdates it (Chang et al 2014).

## MACSE Frameshift Markers Need Post-Processing

MACSE encodes detected frameshifts as `!` (within-codon insertion) and `*` (premature stop) in the nucleotide output. PAML codeml does not accept `!`: codeml 4.10.10 stops with "Error in sequence data file" (and still exits 0), so replace `!` before every codeml run; HyPhy `BUSTED`/`MEME` interpret `!` as `N` but `*` triggers a parse error mid-sequence. Standard cleanup before downstream analysis:

**Goal:** Convert MACSE frameshift and stop markers into formats that PAML and HyPhy will parse correctly.

**Approach:** Replace `!` with gaps for PAML, and use MACSE's own `exportAlignment` sub-program to replace both frameshift codons and internal stops for HyPhy-safe output.

```bash
# Convert frameshift markers to gaps (PAML)
sed -e '/^>/!s/!/-/g' aligned_nt.fasta > aligned_paml.fasta
# Replace frameshift codons and internal stops (HyPhy); checked on MACSE 2.07
java -jar macse_v2.jar -prog exportAlignment -align aligned_nt.fasta \
    -codonForFinalStop --- -codonForInternalStop NNN \
    -codonForInternalFS --- -charForRemainingFS - \
    -out_NT aligned_hyphy.fasta -out_AA aligned_hyphy_aa.fasta
```

`-codonForInternalStop NNN` replaces internal stops with `NNN`; `-codonForInternalFS ---` and `-charForRemainingFS -` remove the `!` frameshift markers, which `exportAlignment` otherwise leaves in place. Without this step, the dN/dS run silently produces results that do not correspond to the alignment shown.

## References

- Chang JM, Di Tommaso P, Notredame C. 2014. TCS: a new multiple sequence alignment reliability measure. MBE 31:1625-1637.
