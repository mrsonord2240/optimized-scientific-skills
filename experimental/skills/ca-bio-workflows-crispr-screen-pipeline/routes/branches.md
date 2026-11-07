# Specialized screen designs

These designs leave this pipeline after counting. Use the named sibling Skill.

| Design | Skill |
|--------|-------|
| Single-cell Perturb-seq, CROP-seq, Multiome (SCEPTRE, Mixscape, Pertpy) | crispr-screens/perturb-seq-analysis |
| Combinatorial paralog, Cas12a Inzolia (interaction scoring) | crispr-screens/combinatorial-screens |
| Base-editor variant-function | crispr-screens/base-editing-analysis, crispr-screens/crispresso-editing |
| Prime-editor variant installation (pegRNA design) | crispr-screens/prime-editing-screens |
| In vivo tumor or immune (per-animal analysis) | crispr-screens/in-vivo-screens |
| Library composition and design | crispr-screens/library-design |

- Single-cell: one guide per cell needs MOI near 0.3; filter or model multi-guide cells.
- In vivo cannot reach 500x coverage: use a focused library (3,000-15,000 guides).
- Hit lists for enrichment or GSEA: pathway-analysis/go-enrichment, pathway-analysis/gsea. MAGeCKFlute gives FluteRRA and FluteMLE dashboards.
