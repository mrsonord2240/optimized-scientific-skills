# When Alignment Is NOT Appropriate (reference)

## When Alignment Is NOT Appropriate

Pick the alignment method by protein identity; below 15% identity DP alignments are statistically indistinguishable from random pairings, so escape to structure or pLM tools.

| Identity (protein) | Recommended approach |
|-------------------|---------------------|
| >= 40% | Any DP aligner; Bio.Align or BLAST is sufficient |
| 25-40% | Use sensitive iterative methods (MMseqs2 iterative, BLASTP with composition-based statistics) |
| 15-25% | Profile-profile (HHsearch, HMMER `phmmer`/`jackhmmer`) |
| < 15% | Structural alignment (Foldseek, TM-align, US-align) or pLM embeddings (TM-Vec, ESM-2 + cosine) -- see structural-alignment |

Length amplifies the signal: 30% identity over 200 residues is far more reliable than 30% over 50.

**When pairwise becomes the wrong tool.** A single DP pairwise alignment is correct for two sequences. For one query against thousands to millions of targets (genome-scale homology search) or for many-vs-many all-by-all (clustering, ortholog detection), the right tool is profile- or k-mer-indexed search, not iterated DP:
- BLASTP / DIAMOND -- standard query-vs-database baseline
- MMseqs2 (Steinegger & Soding 2017 Nat Biotech) -- ~400x faster than PSI-BLAST at higher sensitivity; iterative profile mode `--num-iterations 3` matches PSI-BLAST and approaches `jackhmmer`
- MMseqs2-GPU (Kallenborn et al 2025 Nat Methods 22:2024; Mirdita co-author) -- GPU-accelerated; ~177x faster than `jackhmmer` for single queries on one NVIDIA L40S; use when GPU is available and the dataset is sensitivity-bound
- jackhmmer (HMMER) -- gold standard for distant homology when run to convergence; slow but the most sensitive non-structural method
- Foldseek -- escape to structural search when both query and database have predicted structures (see `alignment/structural-alignment`)

```bash
# one query vs a target FASTA (or MMseqs2 DB); --num-iterations 3 = profile-iterated, PSI-BLAST-like (MMseqs2 18.8cc5c)
mmseqs easy-search query.fa targets.fa hits.tsv tmp --num-iterations 3 -s 7.5 --format-output "query,target,pident,alnlen,evalue,bits"
jackhmmer -N 3 --tblout hits.tbl query.fa targets.fa                # HMMER 3.4; hit table columns: target, ..., full-sequence E-value, score
hhsearch -i query.a3m -d /path/to/hhsuite_db -o query.hhr            # HH-suite 3.3.0 profile-profile; needs a downloaded HH-suite database (below)
```

HBA_HUMAN against 8 UniProt globins: MMseqs2 ranks HBB_HUMAN at 115 bits, E 3e-34; jackhmmer returns all 8 with E < 1e-63. `hhsearch` itself was run (not just checked against `--help`) against a database built from the same 8 sequences (`hhmake` per sequence, `ffindex_build` for the `_hhm`/`_a3m` files, `cstranslate -f` for `_cs219`): querying with HBA_HUMAN ranks itself Prob 100.0/E 1.8e-93, the other alpha globins 100.0 (E 4.6e-77, 6.4e-46), beta-globin orthologs 99.7 (E ~1e-25), myoglobins 96-97 (E ~1e-9) -- probability tracks phylogenetic distance as expected. A real search needs a pre-built database (e.g. `pdb70`, `uniclust30`) fetched with HH-suite's own `hhsuitedb.py`/download scripts, not this toy one.

Other failure modes:
- **Non-homologous sequences**: All DP aligners return an alignment regardless of homology. E-value or bit score is the homology gate, not the existence of an alignment.
- **Repetitive sequences**: Tandem repeats produce ambiguous, artifactually high-scoring alignments; mask first.

## References

- Rost B. 1999. Twilight zone of protein sequence alignments. Prot Eng 12:85-94.
- Steinegger M, Soding J. 2017. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotech 35:1026-1028.
- Kallenborn F, Chacon A, Hundt C, Sirelkhatim H, Didi K, Cha S, Dallago C, Mirdita M, Schmidt B, Steinegger M. 2025. GPU-accelerated homology search with MMseqs2. Nat Methods 22(10):2024-2027.
