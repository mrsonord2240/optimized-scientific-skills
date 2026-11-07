# GO enrichment problems

| Symptom | Cause | Fix |
|---------|-------|-----|
| `No gene can be mapped` | `keytype` or OrgDb does not match the IDs | Match `keytype=` to `keytypes(OrgDb)`; check `orgdb=` species |
| Empty table | `pvalueCutoff` applies to `p.adjust`; universe too large; IDs lost | Rerun with `pvalueCutoff=1 qvalueCutoff=1` to inspect, fix the universe, read the conversion rate in `results.csv.log` |
| Broad terms such as "cellular process" on top | Ranked by p alone; large terms have power at tiny fold change | Read `FoldEnrichment`; keep `minGSSize=10`, `maxGSSize=500` |
| Many near-identical ancestor terms | GO true-path propagation | Keep `simplify=true`, per ontology |
| Terms for tissue-restricted or lowly expressed genes | Whole-genome or default universe | Give the tested-gene universe |
| Long-gene terms (ECM, adhesion) up, ribosomal down in RNA-seq | Detection-power bias | `routes/length-bias.md` |
| Count looks too large | One-to-many ID maps | The scripts deduplicate; check your own pre-mapping |
| Description column shows IDs | Not readable | The script sets `readable=TRUE`; for your own object use `setReadable(ego, OrgDb, 'ENTREZID')` |
| Result is MF when BP was expected | `enrichGO` default is MF | Pass `ont=` explicitly |

Done when the symptom's fix has been applied and the script rerun.
