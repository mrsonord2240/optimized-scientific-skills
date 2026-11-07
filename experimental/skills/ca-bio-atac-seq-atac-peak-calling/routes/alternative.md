# Joint or HMM-based calling

```bash
macs3 hmmratac -i sample.bam -f BAMPE -n sample_hmm --outdir hmm/
```

Use HMM calling only with strong fragment-size periodicity and sufficient nuclear depth; otherwise use standard MACS calling. For Genrich joint-replicate mode use name-sorted BAMs and match the mitochondrial contig name; don't deduplicate twice. Commands and caller tradeoffs are in `references/callers.md`. For single-cell pseudobulk use `bio-atac-seq-single-cell-atac`.

Done when the requested peak file and method-specific quality checks are reported.
