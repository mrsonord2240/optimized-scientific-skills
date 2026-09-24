# SOFT/MINiML formats and legacy resources

## SOFT vs MINiML

| Format | Content | Parser support |
|---|---|---|
| **SOFT** (`*_family.soft.gz`) | Plain-text, key=value style | GEOparse (Python), GEOquery (R), Entrez Direct |
| **MINiML** (`*_family.xml.tgz`) | XML-structured | GEOparse, GEOquery, custom XML |

Both contain the same content. SOFT is the legacy, MINiML the XML successor. GEOparse handles SOFT well; for very large series (1000+ samples) MINiML's XML structure is slower to parse.

## GEOmetadb status

GEOmetadb (Zhu 2008) was a SQLite mirror of GEO metadata enabling fast SQL queries. **Unmaintained since 2020**; downloads still work but data is stale. Modern replacement: pysradb (`pysradb gse_to_srp`, `pysradb metadata`) covers most of the GEO->SRA mapping; for full GEO queries fall back to Entrez gds.

## ArrayExpress -> BioStudies migration (2020)

ArrayExpress (EMBL-EBI's microarray archive, mirroring GEO) was migrated into BioStudies in 2020. Old `E-MTAB-####` accessions still resolve but the API moved:

| Old (pre-2020) | New (BioStudies) |
|---|---|
| `https://www.ebi.ac.uk/arrayexpress/...` | `https://www.ebi.ac.uk/biostudies/...` |
| ArrayExpress REST | BioStudies REST: `https://www.ebi.ac.uk/biostudies/api/v1/...` |

For new workflows, use BioStudies. For legacy ArrayExpress URLs in old papers, redirect via BioStudies.
