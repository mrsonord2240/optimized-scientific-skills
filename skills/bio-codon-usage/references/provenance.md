# Provenance and references

## Provider identity

- Provider: GPTomics
- Repository: `GPTomics/bioSkills`
- Commit: `d91ed3d563019e649dc854c56ccd62551359488a`
- Source path: `sequence-manipulation/codon-usage`
- Source subtree: `3874a00444ef0856451fece587c051cfc34aaf80`
- License: MIT; see [LICENSE](../LICENSE)

The provider's `SKILL.md` and `usage-guide.md` were consolidated into a concise
entrypoint plus routed references. Normalization moved the three reusable
examples from `examples/` to `scripts/` with the exact source blobs below.
Readiness repair then modified those normalized copies to share executable CDS
validation and safe CAI wrappers; the upstream blobs remain the provenance
anchors, not claims that the current files are byte-identical.

| Normalized path | Upstream path | Upstream Git blob |
|---|---|---|
| `scripts/basic_analysis.py` | `examples/basic_analysis.py` | `9bd59859ec36d1c9603d25a8bc58be99ec1d5285` |
| `scripts/cai_optimization.py` | `examples/cai_optimization.py` | `fe608603f49f73fccf3f4844887632af0cec16ce` |
| `scripts/rscu_analysis.py` | `examples/rscu_analysis.py` | `4cc406ef426567b8554c0f7f51c11bbf1b8ebd5c` |

`scripts/codon_utils.py` and `tests/test_codon_usage.py` are optimization-layer
additions. The provider's simplified reciprocal-homozygosity helper was
removed from the advertised surface after bounded comparisons showed it was
not Wright/codonW-compatible Nc.

## Scientific references retained from the provider

- Sharp PM, Li WH. 1987. The codon adaptation index: a measure of
  directional synonymous codon usage bias and its potential applications.
  *Nucleic Acids Research* 15(3):1281-1295.
- dos Reis M, Savva R, Wernisch L. 2004. Solving the riddle of codon usage
  preferences: a test for translational selection. *Nucleic Acids Research*
  32(17):5036-5044.
- Tuller T, Carmi A, Vestsigian K, et al. 2010. An evolutionarily conserved
  mechanism for controlling the efficiency of protein translation. *Cell*
  141(2):344-354.
- Wright F. 1990. The effective number of codons used in a gene. *Gene*
  87(1):23-29.

## Related skill routes retained from the provider

- `transcription-translation` for translation and genetic-code selection.
- `sequence-properties` for GC and sequence-composition analysis.
- `sequence-io/read-sequences` for parsing reference CDS from FASTA or
  GenBank.
- `database-access/entrez-fetch` for retrieving candidate reference genes.

