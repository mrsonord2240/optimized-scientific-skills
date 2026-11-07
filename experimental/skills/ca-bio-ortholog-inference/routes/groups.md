# OrthoDB and OMA

```bash
python scripts/query.py orthodb-search --query "tumor protein p53" --taxid 9606 --out groups.json
```

OrthoDB search is full text, not a gene lookup; inspect names and verify group membership before interpreting target members. Run `python scripts/query.py orthodb-members --group GROUP_ID --taxid 10090 --out members.json` after choosing a group. For OMA use `python scripts/query.py oma --protein P04637 --rel-type 1:1 --out oma.json`; a 404 or empty result is represented as no call by the source helper, whereas exhausted request failures raise. The client also exposes group membership checks, HOG lookup and HOG members for input-specific use; inspect `scripts/ortholog_clients.py`. Resolve identifiers consistently, e.g. using `bio-uniprot-access`.

Done when a verified group or protein lookup has returned interpretable members, or an empty result is explicitly recorded.
