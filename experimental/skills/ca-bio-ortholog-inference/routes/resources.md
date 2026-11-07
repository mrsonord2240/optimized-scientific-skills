# KEGG, PANTHER and eggNOG

```bash
python scripts/query.py kegg --species hsa --gene 7157 --out ko.json
```

The output maps a KEGG gene to KO IDs and lists their member genes; the source example `examples/kegg_orthology.py` also retrieves KO descriptions. KEGG responses are text, and its usage terms vary by use; verify current terms. For PANTHER, the source gives `curl -G 'https://pantherdb.org/services/oai/pantherdb/ortholog/matchortho' --data-urlencode 'geneInputList=P04637' --data-urlencode 'organism=9606' --data-urlencode 'targetOrganism=10090' --data-urlencode 'orthologType=LDO' -o panther.json`; replace IDs after checking the service contract. The source found the eggNOG scripted API unusable; no working API workflow is supplied. Use eggNOG-mapper for local batch annotation or its web interface and adapt to the installed version. HomoloGene is retired and is not a live route.

Done when requested mappings are inspected with resource, release and any missing or failed calls recorded.
