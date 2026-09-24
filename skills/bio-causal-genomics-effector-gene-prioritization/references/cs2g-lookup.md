Moved out of `SKILL.md`: read when looking up cS2G SNP-to-gene scores for credible-set variants.

## cS2G Lookup

**Goal:** Read heritability-calibrated SNP-to-gene scores for a locus's credible-set variants from the pre-computed cS2G table; no install.

**Approach:** Download `cS2G_1000GEUR.zip` (95 MB, hg19, rsID-keyed; per-chromosome `cS2G.<chr>.SGscore.gz` with columns SNP, GENE, cS2G, INFO) from zenodo.org/records/7754032, then run `examples/cs2g_lookup.py` on the credible-set rsIDs. It prints per-SNP gene scores plus a per-gene sum across the queried SNPs (checked 2026-09-21). Lift over first if the credible set is hg38.

```bash
curl -L -o cS2G_1000GEUR.zip https://zenodo.org/api/records/7754032/files/cS2G_1000GEUR.zip/content
python examples/cs2g_lookup.py cS2G_1000GEUR.zip 1 rs11206509 rs10788994   # chr, then rsIDs
```

cS2G scores for one SNP sum to at most 1 across its genes; `INFO` names the constituent strategies that fired (Promoter, ABC, EpiMap, Roadmap, GTeX_Finemapped, eQTLGen_Finemapped). Use it as one evidence stream alongside L2G, not as a substitute (see Reconciliation).
