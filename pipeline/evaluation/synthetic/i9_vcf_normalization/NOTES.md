# i9_vcf_normalization (NOT agent-visible)

## Pitfall and source
Callsets from different callers write the same variant differently: indels not left-aligned, over-padded with extra reference bases, several ALTs in one record. Comparing raw CHROM/POS/REF/ALT undercounts concordance. Fix: split multi-allelics, trim, left-align against the reference, then compare (Tan et al. 2015 Bioinformatics 31:2202, `vt normalize`; `bcftools norm -m-any -f`). Design came from the brief only, no Skill text was read.

## Construction (generate.py, seed 20260919, deterministic, verified by md5 across two runs)
- Reference: ctg1 42 kb, ctg2 31 kb, ctg3 25 kb, random sequence with an embedded repeat tract every ~200 bp (35% homopolymer 6-14, 30% dinucleotide x5-12, 20% trinucleotide x4-8, 15% tetranucleotide x3-7), flanks forced non-periodic so the tract extent is unambiguous. `.fai` verified identical to `samtools faidx`.
- Pool of 297 distinct variants (by left-aligned minimal form): 204 both (C), 52 A-only, 41 B-only. ~58% SNV, 7 MNV (2-3 bp, first and last bases differ so they never trim to SNVs), rest ins/del; 75% of singleton indels sit inside repeat tracts.
- Caller B: canonical left-aligned minimal record per variant, biallelic (245 records).
- Caller A (239 records): each indel shifted to its rightmost equivalent position; about 45% of indel records and 50% of MNV/multi records additionally over-padded with 0-5 left / 0-4 right reference bases. 16 multi-allelic records (SNV+SNV, SNV+SNV+SNV, ins+del of one unit at a tract start, SNV + indel sharing the anchor base). Sites are combined over a union window, so alleles are padded/rewritten.
- Multi-allelic site mix: 10 sites both alts shared, 4 sites one alt shared + one A-only, 2 sites where B has both alts and A only one (A writes biallelic, B two records at one POS), 1 site with both alts A-only, 1 triallelic (2 shared + 1 A-only). So whole-record matching is wrong in several distinct ways.
- Decoys (15 pairs): an A-only and a B-only variant with identical raw CHROM:POS in the two files but different variants (6 same-base SNVs with different ALT; 5 A right-aligned unpadded indel vs B SNV on its anchor base; 4 B left-aligned indel vs A SNV on its anchor base). A raw identical CHROM:POS:REF:ALT record cannot be a false concordance (identical strings normalise identically), so decoys are position-level.
- Truth is the set assignment; generate.py asserts that its own Tan-style normalisation of the written files reproduces (204, 52, 41).

## Verification
| method | both / A-only / B-only |
| --- | --- |
| constructed truth | 204 / 52 / 41 |
| solvers.correct (independent trim + left-align, per allele) | 204 / 52 / 41 |
| `bcftools norm -m-any -f` + `bcftools isec -c none` (bioeval-env, bcftools 1.13) | 204 / 52 / 41 |
| trap raw record compare | 104 / 135 / 141 |
| trap CHROM:POS only | 127 / 112 / 105 |
| trap split, no trim/align | 107 / 149 / 138 |
| trap split+trim, no left-align | 138 / 118 / 107 |
| trap normalise, no split (record level; matches real `bcftools norm -f` w/o `-m` + isec: 178 / 61 / 67) | 178 / 61 / 67 |
| trap first ALT only | 193 / 46 / 52 |

Smallest trap deviation is 6 (first-ALT-only, A-only column) and 11 on n_both; every trap fails on at least two fields.

## Tolerance
Exact (absolute 0). The answer is an integer count fully determined by the stated convention; two independent implementations (Python, bcftools) agree exactly, and there is no method variant (vt, bcftools, GATK LeftAlignAndTrimVariants, hand-rolled) that should legitimately differ. Any non-zero tolerance would let the smallest trap (first ALT only, off by 11/6/11) close in.

## Task wording deviation from the brief
The brief said to state "identical normalised representation (left-aligned, minimal, one alt per record)". That names the method, which the contract forbids, so the task instead defines identity semantically ("same change to the reference sequence at the same place, regardless of how each caller wrote it") and states one-variant-per-ALT, ignore GT/QUAL/DP, no MNV decomposition, distinct-variant counting. These are equivalent for this data because every pool variant is a distinct left-aligned form.

## Known limits / ambiguity / gaming
- Counting a multi-allelic record per record rather than per ALT is explicitly ruled out in the text, but a solver skimming may still do it (that is trap_norm_no_split territory).
- An agent that reasons at the haplotype level (apply each ALT to the reference and compare sequences) gets the same answer; that is legitimate.
- MNVs: an agent that decomposes them would change counts; the text forbids it. Only 7 MNVs, so a decomposition mistake alone changes counts by a small amount but still fails at tolerance 0 (wanted: exact answer).
- Gameable only by knowing the generator; no file in data/ marks truth. Sample ID `SAMPLE1` and `##source=callerA/B` are the only labels.
- bcftools 1.13 (sandbox) tested; newer versions behave identically for this input (no `-m` corner cases: alleles are simple).
- Not a benchmark score; fabricated data.
