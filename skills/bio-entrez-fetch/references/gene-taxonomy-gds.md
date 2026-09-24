# Gene, taxonomy and GEO DataSets (gds) retrieval

Read when fetching from `gene`, `taxonomy` (lineage by TXID) or `gds`.

Assumes the Required Setup block in `SKILL.md` (`Entrez.email`, imports, `expect_start`).

### gene

| rettype | retmode | Returns | Use when |
|---|---|---|---|
| `gene_table` | `text` | Tabular per-transcript layout | Exon coordinates |
| `xml` | `xml` | Full Entrez Gene XML | Everything else — name, synonyms, GeneRIFs, locus |

### taxonomy

| rettype | retmode | Returns | Use when |
|---|---|---|---|
| `xml` | `xml` (default) | TaxNode XML | Lineage, parent, common name |

### gds (GEO)

| rettype | retmode | Returns | Use when |
|---|---|---|---|
| (default — no rettype) | `text` | Plaintext SOFT-style summary | Quick metadata; for full series matrix go to FTP |

EFetch for GDS records is intentionally minimal — full GEO downloads go via the FTP mirror or `GEOparse`. See `geo-data` skill.

### Taxonomy lineage by TXID

```python
def lineage(txid):
    h = Entrez.efetch(db='taxonomy', id=str(txid), retmode='xml')
    record = Entrez.read(h)[0]; h.close()
    return record['Lineage'], record['ScientificName']
```

Given a species name instead of a TXID, ESearch the `taxonomy` db first to resolve it — never assume a TXID from memory.
