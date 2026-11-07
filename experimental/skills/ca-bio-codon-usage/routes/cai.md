# CAI and max-CAI example

```bash
python scripts/cai_optimization.py
```

The source example uses toy reference and query sequences; it does not accept files or provide a validated host reference. For meaningful CAI, adapt it to validated in-frame CDS and an expression-biased reference from the relevant host, using the same genetic code for the index and guarded score. Record reference provenance and excluded codons. Biopython behavior is version-sensitive; inspect [compatibility](../references/compatibility-and-errors.md). The example also produces a naive single-codon max-CAI draft through the source table-aware helper; a high score does not establish expression or design readiness. Its omitted constraints and tie behavior are described in [constraints](../references/optimization-constraints.md). tAI and a standard Nc estimator are not implemented.

**Done when:** the reference, table, validation, score and limitations are recorded; any example output is clearly identified as an unscreened draft.
